// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "HyPhyStudio",
    platforms: [.macOS(.v14)],
    products: [
        .executable(name: "HyPhyStudio", targets: ["HyPhyStudio"]),
    ],
    targets: [
        .target(
            name: "HyPhyBridge",
            path: "Sources/HyPhyBridge",
            publicHeadersPath: "include",
            cxxSettings: [
                .unsafeFlags([
                    "-std=c++20",
                    "-O3",
                    "-I../../include",
                    "-I../../build/_deps/eigen3-src",
                    "-I../../_deps/eigen3-src"
                ])
            ]
        ),
        .executableTarget(
            name: "HyPhyStudio",
            dependencies: ["HyPhyBridge"],
            path: "Sources/HyPhyStudio",
            resources: [
                .process("Resources")
            ],
            swiftSettings: [
                .interoperabilityMode(.Cxx)
            ]
        )
    ]
)
