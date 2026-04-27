# 🧬 Plasmid Shredder & Cocktail Optimizer

**Plasmid Shredder** is an interactive bioinformatics tool designed to automate the search for ideal restriction enzyme cocktails. It helps molecular biologists identify combinations of enzymes (up to 3) that will digest a circular plasmid into fragments within a specific size range, while simultaneously optimizing for buffer compatibility and experimental conditions.

---

## 🚀 Live Demo
You can view and use the app here: **[INSERT YOUR STREAMLIT CLOUD LINK HERE]**

---

## ✨ Key Features

* **Combinatorial Digestion Search:** Uses a highly optimized algorithm to test millions of enzyme combinations in seconds.
* **NEB Compatibility Scoring:** Automatically calculates the best shared buffer (r1.1, r2.1, r3.1, or CutSmart) and scores cocktails based on:
    * **Activity Levels:** Penalizes combinations where an enzyme has <100% activity.
    * **Temperature Matching:** Ensures all enzymes in the cocktail operate at the same incubation temperature (e.g., 37°C).
    * **Star Activity Alerts:** Identifies and penalizes buffers that risk non-specific cutting.
* **Fragment Filtering:** Set strict minimum and maximum fragment size boundaries to ensure your "shredded" DNA is perfect for downstream applications like NGS, Gibson Assembly, or capillary electrophoresis.
* **Interactive Dashboard:** Built with Streamlit for a clean, easy-to-use web interface.

## 🛠️ Installation & Local Usage

If you wish to run this tool locally on your own machine:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/](https://github.com/)[YOUR-USERNAME]/plasmid-shredder.git
   cd plasmid-shredder
