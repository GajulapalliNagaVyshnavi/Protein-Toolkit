from setuptools import setup, find_packages

# Read the README file for long description
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "A comprehensive toolkit for protein sequence analysis"

setup(
    name="Protein_Toolkit",
    version="1.0.0",
    author="JoVy",
    author_email="nvyshnavi36@gmail.com",
    description="A comprehensive toolkit for protein sequence analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/GajulapalliNagaVyshnavi/Protein-Toolkit",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",  # Changed from "==3.10" to be more flexible
    install_requires=[
        # Core dependencies
        "pandas>=2.0.0",
        "numpy>=2.0.0",
        "biopython>=1.80",
    ],
    extras_require={
        # For development
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
        ],
        # For deep learning features (optional)
        "ml": [
            "torch>=2.9.0",
            "transformers>=4.57.0",
            "fair-esm>=2.0.0",
        ],
        # All dependencies
        "all": [
            "torch>=2.9.0",
            "transformers>=4.57.0",
            "fair-esm>=2.0.0",
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
        ]
    },
)