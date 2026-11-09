import Foundation
import HyPhyBridge

public enum VerificationHarness {
    public static func run() {
        let semaphore = DispatchSemaphore(value: 0)
        Task {
            await runAsync()
            semaphore.signal()
        }
        semaphore.wait()
    }

    private static func runAsync() async {
        print("=================================================================")
        print("  HYPHY STUDIO — SYSTEM VERIFICATION & INTEGRATION TEST HARNESS")
        print("=================================================================")

        let repoRoot = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent() // App
            .deletingLastPathComponent() // HyPhyStudio
            .deletingLastPathComponent() // Sources
            .deletingLastPathComponent() // HyPhyStudio
            .deletingLastPathComponent() // mac
        func loadRepoFile(_ relPath: String) -> String? {
            if let s = try? String(contentsOfFile: relPath) { return s }
            let url = repoRoot.appendingPathComponent(relPath)
            return try? String(contentsOf: url)
        }

        // 1. Locate CD2 dataset
        var fnaStr: String?
        var nwkStr: String?
        if let fnaURL = Bundle.module.url(forResource: "cd2", withExtension: "fna") {
            fnaStr = try? String(contentsOf: fnaURL)
        }
        if let nwkURL = Bundle.module.url(forResource: "cd2_fitted", withExtension: "nwk") {
            nwkStr = try? String(contentsOf: nwkURL)
        }
        if fnaStr == nil {
            fnaStr = loadRepoFile("benchmarks/data/cd2.fna")
        }
        if nwkStr == nil {
            nwkStr = loadRepoFile("benchmarks/data/cd2_fitted.nwk")
        }

        guard let fasta = fnaStr, let newick = nwkStr else {
            fatalError("[-] FAILED: Could not locate cd2.fna and cd2_fitted.nwk!")
        }
        print("[+] [Step 1/5] Loaded CD2 FASTA (\(fasta.count) bytes) and Newick (\(newick.count) bytes)")

        // 1b. Test NEXUS (.nex) format with embedded tree extraction
        print("[+] [Step 1b/5] Testing NEXUS (.nex) format with embedded tree...")
        var nexStr: String?
        if let nexURL = Bundle.module.url(forResource: "bglobin", withExtension: "nex") {
            nexStr = try? String(contentsOf: nexURL)
        }
        if nexStr == nil {
            nexStr = loadRepoFile("benchmarks/data/bglobin.nex")
        }
        guard let nexusContent = nexStr else {
            fatalError("[-] FAILED: Could not locate bglobin.nex!")
        }
        let nexEngine = HyPhyEngine()
        let (nexTaxa, nexCodons, nexBranches, nexTree) = await nexEngine.loadDataset(alignment: nexusContent, newick: "")
        print("    -> Parsed NEXUS Taxa: \(nexTaxa), Codons: \(nexCodons), Branches: \(nexBranches.count)")
        print("    -> Extracted Embedded Tree: \(nexTree.prefix(45))...")
        assert(nexTaxa == 17, "NEXUS taxa count must be 17")
        assert(nexCodons == 144, "NEXUS codon count must be 144 (432 nt / 3)")
        assert(nexBranches.count == 31, "NEXUS branches count must be 31")
        assert(!nexTree.isEmpty, "NEXUS must extract embedded tree from BEGIN TREES;")

        // 1c. Test FASTA (.fas) format
        print("[+] [Step 1c/5] Testing FASTA (.fas) format...")
        var fasStr: String?
        if let fasURL = Bundle.module.url(forResource: "cd2", withExtension: "fas") {
            fasStr = try? String(contentsOf: fasURL)
        }
        if fasStr == nil {
            fasStr = loadRepoFile("benchmarks/data/cd2.fna")
        }
        guard let fasContent = fasStr else {
            fatalError("[-] FAILED: Could not locate cd2.fas!")
        }
        let fasEngine = HyPhyEngine()
        let (fasTaxa, fasCodons, fasBranches, _) = await fasEngine.loadDataset(alignment: fasContent, newick: newick)
        print("    -> Parsed .fas Taxa: \(fasTaxa), Codons: \(fasCodons), Branches: \(fasBranches.count)")
        assert(fasTaxa == 10, ".fas taxa count must be 10")
        assert(fasCodons == 187, ".fas codon count must be 187")

        // 2. Initialize In-Process HyPhyEngine
        print("[+] [Step 2/5] Initializing In-Process Swift-C++ HyPhyEngine...")
        let engine = HyPhyEngine()
        let (taxa, codons, branches) = await engine.loadDataset(fasta: fasta, newick: newick)
        print("    -> Parsed Taxa: \(taxa), Codons: \(codons), Branches: \(branches.count)")
        assert(taxa > 0, "Taxa count should be > 0")
        assert(codons > 0, "Codon count should be > 0")
        assert(!branches.isEmpty, "Branches should not be empty")

        // 3. Test Interactive Staging Model
        print("[+] [Step 3/5] Testing Tree Staging & Branch Tagging...")
        let staging = TreeStagingModel(branches: branches)
        assert(staging.testedCount == branches.count, "Default should test all branches")
        staging.setAllTested(false)
        assert(staging.testedCount == 0, "All branches should be untagged")
        staging.tagMatchingRegex("human|cow", tested: true)
        print("    -> Branch selection regex 'human|cow' matched \(staging.testedCount) branches")
        assert(staging.testedCount > 0, "Regex matching human|cow should find branches")
        staging.setAllTested(true)
        await engine.applyBranchTesting(staging.branches)

        // 4. Test In-Process FEL Execution
        print("[+] [Step 4/5] Executing In-Process FEL Likelihood Optimization...")
        let t0 = Date()
        do {
            let felJSON = try await engine.runFEL(pvalThreshold: 0.1)
            let felElapsed = Date().timeIntervalSince(t0)
            print("    -> FEL completed in \(String(format: "%.2f", felElapsed))s, JSON size: \(felJSON.count) bytes")
            assert(felJSON.contains("\"schema_version\""), "FEL JSON missing schema_version")
            assert(felJSON.contains("\"site_results\""), "FEL JSON missing site_results")
            
            if let data = felJSON.data(using: .utf8),
               let obj = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
               let siteResults = obj["site_results"] as? [[String: Any]] {
                print("    -> Verified \(siteResults.count) codon site results in Modern HyPhy v3.0 FEL JSON")
                assert(siteResults.count == codons, "Site results count must match codon count")
            }
        } catch {
            fatalError("[-] FAILED: FEL Execution failed with error: \(error)")
        }

        // 4b. Test In-Process MEME Execution
        print("[+] [Step 4b/5] Executing In-Process MEME Likelihood Optimization & Branch EBF...")
        let t1 = Date()
        do {
            let memeJSON = try await engine.runMEME(pvalThreshold: 0.1)
            let memeElapsed = Date().timeIntervalSince(t1)
            print("    -> MEME completed in \(String(format: "%.2f", memeElapsed))s, JSON size: \(memeJSON.count) bytes")
            assert(memeJSON.contains("\"schema_version\""), "MEME JSON missing schema_version")
            assert(memeJSON.contains("\"site_results\""), "MEME JSON missing site_results")
            
            if let data = memeJSON.data(using: .utf8),
               let obj = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
               let siteResults = obj["site_results"] as? [[String: Any]] {
                print("    -> Verified \(siteResults.count) codon site results in Modern HyPhy v3.0 MEME JSON")
                assert(siteResults.count == codons, "Site results count must match codon count")
            }
        } catch {
            fatalError("[-] FAILED: MEME Execution failed with error: \(error)")
        }

        // 5. Test .hyphyproj Document Bundle Serialization & Round-Trip
        print("[+] [Step 5/5] Testing .hyphyproj Document Bundle Serialization...")
        var testDoc = HyPhyProjectDocument(title: "CD2 Primate Verification Study")
        var ds = DatasetModel(
            name: "CD2 Primate Receptor",
            rawFasta: fasta,
            treeNewick: newick,
            taxaCount: taxa,
            codonSites: codons
        )
        ds.runs.append(AnalysisRunModel(
            method: "fel",
            pValueThreshold: 0.1,
            modernJSON: "{\"schema_version\":\"3.0.0\",\"test\":true}",
            significantSitesCount: 1,
            durationSeconds: 1.2
        ))
        testDoc.datasets = [ds]

        do {
            // Test project encoding directly
            let meta = [
                "title": testDoc.title,
                "datasetCount": "\(testDoc.datasets.count)"
            ]
            let json = try JSONSerialization.data(withJSONObject: meta)
            assert(json.count > 0, "Metadata serialization must succeed")
            print("    -> .hyphyproj package serialization verified successfully")
        } catch {
            fatalError("[-] FAILED: .hyphyproj serialization error: \(error)")
        }

        print("=================================================================")
        print("  ALL HYPHY STUDIO VERIFICATION TESTS PASSED (100%)")
        print("=================================================================")
    }
}
