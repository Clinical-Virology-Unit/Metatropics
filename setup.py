from setuptools import setup, find_packages
from pathlib import Path

readme_path = Path(__file__).resolve().parent / "README.md"
with readme_path.open("r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="virasign",
    version="0.0.6",
    author="Virasign Team",
    description="Virasign: Viral Read ASSIGNment from nanopore sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(where="nf-metatropics/assets/submission"),
    package_dir={"": "nf-metatropics/assets/submission"},
    python_requires=">=3.7,<3.14",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "virasign=virasign.virasign:main",
            "metatropics_samplesheet=metatropics_samplesheet.samplesheet:main",
            "metatropics-samplesheet=metatropics_samplesheet.samplesheet:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)