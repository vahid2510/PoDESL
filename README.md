# PoDESL — Portable Domain-Specific Language for Finite Element Simulation

> **Portable, Declarative, and Solver-Driven**  
> A compact, high-fidelity DSL and runtime ecosystem for structural, thermal, and multiphysics finite element simulations, designed for rapid prototyping, academic research, and industrial-grade analysis.

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)


---

## 🚀 Project Overview

**PoDESL (Portable Domain-Specific Language)** is not merely a scripting tool; it is a **full-stack computational engineering platform**. Born from a need for rapid, safe, and expressive simulation prototyping, PoDESL has matured into a robust ecosystem encompassing:

1.  **A Declarative Domain-Specific Language (DSL):** Write simulation problems with the clarity of a technical report.
2.  **An Extensive Suite of Solvers:** Over 80 specialized solvers covering structural mechanics, heat transfer, dynamics, and coupled multiphysics.
3.  **A Browser-Based Integrated Development Environment (IDE):** Featuring notebook mode, real-time plotting, VTK export, and Abaqus interoperability.
4.  **A Cross-Platform Packaging Pipeline:** Producing standalone executables for Windows, Linux, and macOS.

This project embodies a philosophy of **"Legibility First, Tooling Parity, Safety Always."** Whether you are an analyst running parametric studies, a student learning FEM, or a developer extending the solver suite, PoDESL provides the tools and documentation to succeed.

---

## 📜 Core Philosophy & Design Principles

PoDESL is built upon three foundational pillars:

### 1. Legibility
A PoDESL file reads like a well-structured technical specification. Problem statements, inputs, solution steps, and reporting directives are grouped into explicit blocks (`GIVEN`, `SOLVE`, `EQUATIONS`, `REPORT`). This structure ensures that even complex simulations remain human-readable and maintainable.

### 2. Composability
DSL snippets can call other DSL files (e.g., study drivers) or export intermediate data that flows into downstream tools like Abaqus, MATLAB, or custom post-processing scripts. This enables the construction of sophisticated workflows from simple, reusable components.

### 3. Tooling Parity
The *exact same* DSL file can be executed via:
*   The command-line interface: `python -m podesl.cli run my_problem.dsl`
*   The browser-based IDE (editor or notebook mode)
*   Automated CI/CD smoke tests
This ensures consistency and eliminates "it works on my machine" issues.

### 4. Safety & Sandboxing
All expressions within the DSL are evaluated in a tightly controlled environment. The Python `eval` function is used with a whitelisted set of globals (NumPy, math, and custom helpers), preventing arbitrary code execution and ensuring the integrity of the simulation environment.

---

## 🧩 Architecture & Components

The PoDESL architecture is modular and well-documented, facilitating easy contribution and extension.

### 1. Parser & Runtime Internals
The core engine parses the DSL text into an Abstract Syntax Tree (AST), transpiles it into executable Python code, and manages its execution.

*   **Parsing Pipeline:** The parser (`podesl/parser.py`) uses a line-reader and block-scanner to identify `PROBLEM` headers and `GIVEN`, `SOLVE`, `EQUATIONS`, and `REPORT` blocks. It enforces consistent indentation and validates syntax.
*   **Transpilation:** The transpiler (`podesl/transpile.py`) converts the AST into a pure Python script. Variables defined in `GIVEN` become assignments, while `REPORT` commands translate into function calls for printing, exporting, or plotting.
*   **Execution Helpers:** The `podesl/execution.py` module handles the safe execution of the generated Python code, capturing stdout/stderr and storing the last solver result for IDE integration.

### 2. Solver Catalogue
PoDESL ships with over 80 solvers, categorized into families:

#### Structural Solvers
*   **Truss Solvers:** `truss2d_linear`, `truss2d_nonlinear`, `truss2d_modal`, `truss3d`.
*   **Frame & Beam Solvers:** `frame2d_static`, `frame2d_pdelta`, `frame2d_transient_newmark`, `frame2d_buckling`, `frame2d_modal`, `frame2d_modal_consistent`, `frame2d_static_nl`, `frame3d_beamcolumn`.
*   **Solid & Shell Solvers:** `elastic3d`, `solid3d_hex8_nonlinear`, `solid3d_hex20`, `shell2d_mindlin_static`, `plate2d_mindlin_modal`, `plate2d_mindlin_static`.

#### Thermal Solvers
*   **1D/2D/3D Conduction:** `heat1d`, `heat2d_q4_static`, `heat2d_steady`, `heat2d_transient`, `heat2d_transient_conv`, `heat3d_h8_transient`, `heat3d_transient`, `heat_axi_transient`.
*   **Transient Analysis:** All thermal solvers support backward Euler, Crank-Nicolson, and explicit schemes with configurable time-stepping.

#### Multiphysics & Study Solvers
*   **Fusion Studies:** `study_fusion` for steady and transient coupling between thermal and mechanical domains.
*   **Monte Carlo & Scenario Studies:** `study_montecarlo` for statistical analysis and `study_scenarios` for comparative benchmarking.
*   **Optimization:** `study_optimize` for gradient-based design optimization.

#### Utility Modules
*   **`bc_utils`:** Boundary condition handling and parsing.
*   **`input_utils`:** Input/output helpers for meshes and results.
*   **`materials`:** Material property management.
*   **`vtk_writer`:** Export of results to VTK format for visualization in ParaView or VisIt.

Each solver is documented in Chapter 8 of the [User Manual](docs/PoDESL User Manual.pdf), detailing its governing equations, input/output structure, and limitations.

### 3. IDE & Workflow Automation
The browser-based IDE provides a rich, interactive experience.

*   **Architecture:** Built with plain JavaScript, it features an Ace editor, Chart.js for plotting, and a Markdown report builder.
*   **Notebook Mode:** Users can switch to a Jupyter-like notebook interface, where each cell can contain a DSL snippet and display its output.
*   **Diagnosis & History:** Run results are recorded server-side, with diagnostics converted into Ace annotations for easy error tracing.
*   **Collaboration & Export:** Results can be copied to the clipboard, exported as `.inp` files for Abaqus, or downloaded as formatted Markdown reports.

### 4. Testing, Build, & Release Process
PoDESL employs rigorous quality assurance practices.

*   **Automated Testing:** Unit tests live under `tests/`, and `tools/smoke_examples.py` executes all example files to ensure no regressions slip in.
*   **Build & Packaging:** The `build_windows.ps1` script creates a PyInstaller one-file executable for Windows. Future releases will include Linux and macOS bundles.
*   **Documentation Pipeline:** The manual is built using LaTeX (`pdflatex` or `latexmk`) and is automatically rendered as part of release candidates.

---

## 🛠️ Getting Started

### Prerequisites
*   Python 3.8 or higher
*   Git

### Installation
Clone the repository and install dependencies:

```bash
git clone https://github.com/vahid2510/PoDESL.git
cd PoDESL
pip install -r requirements.txt
```
### Running Your First Simulation
Via Command Line:
```bash
python -m podesl.cli run examples/truss2d_simple.dsl
```
Via Web IDE:
```bash
python -m podesl.cli --ide
```
Open your browser to http://localhost:5000 and start editing.

### 📚 Documentation
The official reference is the PoDESL User Manual, a comprehensive 200+ page document auto-generated from the source code and examples.

Chapter 1: DSL Overview (Grammar, Syntax, Examples, Style Guide)
Chapter 2: IDE and Workflow Automation (Architecture, Notebook Mode, Export)
Chapter 3: Parser and Runtime Internals (AST, Transpilation, Safety)
Chapter 4: Multiphysics and Study Solvers (Fusion, Monte Carlo, Optimization)
Chapter 5: Structural Solvers (Truss, Frame, Solid, Shell)
Chapter 6: Thermal and Thermo-Mechanical Solvers (Equations, 1D/2D/3D, Coupling)
Chapter 7: Testing, Build, and Release Process
Chapter 8: Detailed Solver Documentation (80+ solvers)
Appendices: Formula Reference, Example Directory
### 🤝 Contributing
Contributions are welcome! Please follow these guidelines:

Read the Manual: Familiarize yourself with the architecture and coding standards outlined in the manual (especially Chapters 1, 3, and 7).
Extend the Parser: To add a new block type or statement, update the parser to recognize the syntax and extend the transpiler to generate the corresponding Python code.
Add a New Solver: Create a new module under podesl/solvers/, implement the physics, and document it thoroughly in Chapter 8 of the manual.
Update Tests: Add unit tests under tests/ and ensure tools/smoke_examples.py passes.
Document Changes: Update the manual, README.md, and any relevant examples.
For major changes, open an issue first to discuss the proposed feature.

### 📊 Use Cases
PoDESL is designed for a wide range of applications:

Academic Research & Teaching: Ideal for demonstrating FEM concepts, conducting benchmark studies, and assigning projects.
Rapid Prototyping: Quickly test new solver algorithms or model configurations before implementing them in larger frameworks.
Industrial Analysis: Perform parametric studies, sensitivity analyses, and preliminary design checks.
Multiphysics Workflows: Model complex interactions between structural, thermal, and fluid domains.
👤 Author & Acknowledgements
PoDESL was conceived, designed, and implemented by Vahid Ahmadi Kharmai.

"This project represents a culmination of years of work in computational mechanics and software engineering. Its goal is to provide a safe, expressive, and powerful tool for engineers and researchers to explore complex physical phenomena."

— Vahid Ahmadi Khormai

📧 Contact: vahid.ahmadi.kharmai@Gmail.com 


📜 License
This project is licensed under the MIT License. See the LICENSE file for details.

© 2025 Vahid Ahmadi Khormai. All rights reserved.
```

