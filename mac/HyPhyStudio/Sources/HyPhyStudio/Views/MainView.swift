import SwiftUI
import UniformTypeIdentifiers

public struct MainView: View {
    @Binding public var document: HyPhyProjectDocument

    @State private var selectedDatasetID: UUID?
    @State private var selectedRunID: UUID?
    @State private var showingRunSheet = false
    @State private var treeStagingModel = TreeStagingModel()
    @State private var stagingMode: StagingTab = .treeStaging

    public enum StagingTab: String, CaseIterable {
        case treeStaging = "Phylogeny Staging"
        case alignmentQC = "Alignment QC"
    }

    public init(document: Binding<HyPhyProjectDocument>) {
        self._document = document
    }

    private var selectedDataset: DatasetModel? {
        document.datasets.first(where: { $0.id == selectedDatasetID })
    }

    private var selectedRun: AnalysisRunModel? {
        selectedDataset?.runs.first(where: { $0.id == selectedRunID })
    }

    public var body: some View {
        NavigationSplitView {
            SidebarView(
                datasets: $document.datasets,
                selectedDatasetID: $selectedDatasetID,
                selectedRunID: $selectedRunID,
                onAddDataset: promptAddDataset,
                onImportSample: loadCD2Sample
            )
            .navigationSplitViewColumnWidth(min: 210, ideal: 240, max: 320)
        } content: {
            if let run = selectedRun {
                // 1. Render Results in Dieter Rams Visualizer (WKWebView)
                VisualizerWebView(resultJSON: run.modernJSON)
                    .ignoresSafeArea(.all, edges: .bottom)
            } else if let dataset = selectedDataset {
                // 2. Pre-flight Staging View (Tree branch tagging + QC)
                VStack(spacing: 0) {
                    Picker("", selection: $stagingMode) {
                        ForEach(StagingTab.allCases, id: \.self) { tab in
                            Text(tab.rawValue).tag(tab)
                        }
                    }
                    .pickerStyle(.segmented)
                    .padding(8)

                    Divider()

                    if stagingMode == .treeStaging {
                        TreeStagingView(model: treeStagingModel)
                    } else {
                        // Alignment QC View
                        VStack(alignment: .leading, spacing: 12) {
                            Text("Alignment Quality Control")
                                .font(.headline)
                            Text("Taxa: \(dataset.taxaCount) • Codon Sites: \(dataset.codonSites) • Genetic Code: \(dataset.geneticCode)")
                                .font(.subheadline)
                                .foregroundStyle(.secondary)

                            ScrollView {
                                Text(dataset.rawFasta.prefix(2000))
                                    .font(.system(size: 11, design: .monospaced))
                                    .frame(maxWidth: .infinity, alignment: .leading)
                                    .padding(10)
                                    .background(Color(nsColor: .controlBackgroundColor))
                                    .clipShape(RoundedRectangle(cornerRadius: 6))
                            }
                        }
                        .padding(16)
                    }
                }
            } else {
                ContentUnavailableView(
                    "No Dataset Selected",
                    systemImage: "dna",
                    description: Text("Select a dataset from the sidebar or click '+' to import an alignment and tree.")
                )
            }
        } detail: {
            InspectorView(dataset: selectedDataset, run: selectedRun)
        }
        .toolbar {
            ToolbarItemGroup(placement: .primaryAction) {
                Button(action: { showingRunSheet = true }) {
                    Label("Run Analysis", systemImage: "play.fill")
                }
                .disabled(selectedDataset == nil)
                .keyboardShortcut("r", modifiers: .command)
            }
        }
        .sheet(isPresented: $showingRunSheet) {
            if let ds = selectedDataset {
                RunSheetView(dataset: ds, treeModel: treeStagingModel) { newRun in
                    if let idx = document.datasets.firstIndex(where: { $0.id == ds.id }) {
                        document.datasets[idx].runs.append(newRun)
                        selectedRunID = newRun.id
                    }
                }
            }
        }
        .onChange(of: selectedDatasetID) { _, newID in
            if let id = newID, let ds = document.datasets.first(where: { $0.id == id }) {
                loadTreeIntoStaging(ds)
            }
        }
        .onAppear {
            if document.datasets.isEmpty {
                loadCD2Sample()
            } else if selectedDatasetID == nil {
                selectedDatasetID = document.datasets.first?.id
            }
        }
    }

    private func loadTreeIntoStaging(_ ds: DatasetModel) {
        let engine = HyPhyEngine()
        Task {
            let (taxa, codons, branches, effectiveTree) = await engine.loadDataset(alignment: ds.rawFasta, newick: ds.treeNewick)
            await MainActor.run {
                if !branches.isEmpty {
                    treeStagingModel = TreeStagingModel(branches: branches)
                }
                if ds.treeNewick.isEmpty && !effectiveTree.isEmpty {
                    if let idx = document.datasets.firstIndex(where: { $0.id == ds.id }) {
                        document.datasets[idx].treeNewick = effectiveTree
                        document.datasets[idx].taxaCount = taxa
                        document.datasets[idx].codonSites = codons
                    }
                }
            }
        }
    }

    private func promptAddDataset() {
        let panel = NSOpenPanel()
        let nexType = UTType(filenameExtension: "nex") ?? .plainText
        let nexusType = UTType(filenameExtension: "nexus") ?? .plainText
        let fasType = UTType(filenameExtension: "fas") ?? .plainText
        let fastaType = UTType(filenameExtension: "fasta") ?? .plainText
        let fnaType = UTType(filenameExtension: "fna") ?? .plainText
        let faType = UTType(filenameExtension: "fa") ?? .plainText
        let nwkType = UTType(filenameExtension: "nwk") ?? .plainText
        let newickType = UTType(filenameExtension: "newick") ?? .plainText
        let treeType = UTType(filenameExtension: "tree") ?? .plainText

        panel.allowedContentTypes = [nexType, nexusType, fasType, fastaType, fnaType, faType, nwkType, newickType, treeType, .plainText]
        panel.allowsOtherFileTypes = true
        panel.allowsMultipleSelection = true
        panel.canChooseFiles = true
        panel.canChooseDirectories = false
        panel.message = "Select Codon Alignment (.nex, .fas, .fasta, .fna) and optional Tree (.nwk)"

        guard panel.runModal() == .OK, !panel.urls.isEmpty else { return }

        var alignmentURL: URL?
        var treeURL: URL?
        var alignmentContent = ""
        var treeContent = ""

        for url in panel.urls {
            let ext = url.pathExtension.lowercased()
            if ext == "nwk" || ext == "newick" || ext == "tree" {
                treeURL = url
                treeContent = (try? String(contentsOf: url)) ?? ""
            } else {
                alignmentURL = url
                alignmentContent = (try? String(contentsOf: url)) ?? ""
            }
        }

        guard let primaryURL = alignmentURL ?? treeURL, !alignmentContent.isEmpty || !treeContent.isEmpty else { return }

        // If only alignment was picked, look in the same directory for a matching tree
        if treeContent.isEmpty, let alnURL = alignmentURL {
            let baseDir = alnURL.deletingLastPathComponent()
            let baseName = alnURL.deletingPathExtension().lastPathComponent
            let candidateNwk = baseDir.appendingPathComponent("\(baseName).nwk")
            let candidateFitted = baseDir.appendingPathComponent("\(baseName)_fitted.nwk")
            if let nwkStr = try? String(contentsOf: candidateNwk) {
                treeContent = nwkStr
            } else if let nwkStr = try? String(contentsOf: candidateFitted) {
                treeContent = nwkStr
            }
        }

        let datasetName = primaryURL.deletingPathExtension().lastPathComponent
        let engine = HyPhyEngine()

        Task {
            let (taxa, codons, branches, effectiveTree) = await engine.loadDataset(alignment: alignmentContent, newick: treeContent)
            await MainActor.run {
                let ds = DatasetModel(
                    name: datasetName,
                    rawFasta: alignmentContent,
                    treeNewick: effectiveTree.isEmpty ? treeContent : effectiveTree,
                    taxaCount: taxa,
                    codonSites: codons
                )
                document.datasets.append(ds)
                selectedDatasetID = ds.id
                if !branches.isEmpty {
                    treeStagingModel = TreeStagingModel(branches: branches)
                }
            }
        }
    }

    private func loadCD2Sample() {
        var fnaStr: String?
        var nwkStr: String?
        var jsonStr: String?

        if let fnaURL = Bundle.module.url(forResource: "cd2", withExtension: "fna") {
            fnaStr = try? String(contentsOf: fnaURL)
        }
        if let nwkURL = Bundle.module.url(forResource: "cd2_fitted", withExtension: "nwk") {
            nwkStr = try? String(contentsOf: nwkURL)
        }
        if let jsonURL = Bundle.module.url(forResource: "cd2_meme", withExtension: "json") {
            jsonStr = try? String(contentsOf: jsonURL)
        }

        // Fallbacks
        let repoRoot = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent() // Views
            .deletingLastPathComponent() // HyPhyStudio
            .deletingLastPathComponent() // Sources
            .deletingLastPathComponent() // HyPhyStudio
            .deletingLastPathComponent() // mac
        if fnaStr == nil {
            fnaStr = (try? String(contentsOfFile: "benchmarks/data/cd2.fna"))
                ?? (try? String(contentsOf: repoRoot.appendingPathComponent("benchmarks/data/cd2.fna")))
        }
        if nwkStr == nil {
            nwkStr = (try? String(contentsOfFile: "benchmarks/data/cd2_fitted.nwk"))
                ?? (try? String(contentsOf: repoRoot.appendingPathComponent("benchmarks/data/cd2_fitted.nwk")))
        }

        guard let fasta = fnaStr, let tree = nwkStr else { return }

        var cd2 = DatasetModel(
            name: "CD2 Primate Receptor",
            rawFasta: fasta,
            treeNewick: tree,
            taxaCount: 10,
            codonSites: 187,
            tags: ["Immune", "Primates", "Benchmark"]
        )

        if let js = jsonStr {
            let run = AnalysisRunModel(
                method: "meme",
                pValueThreshold: 0.10,
                modernJSON: js,
                significantSitesCount: 16,
                durationSeconds: 2.62,
                status: "completed"
            )
            cd2.runs.append(run)
        }

        document.datasets.append(cd2)
        selectedDatasetID = cd2.id
        selectedRunID = cd2.runs.first?.id
    }
}
