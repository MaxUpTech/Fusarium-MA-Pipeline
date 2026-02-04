from setuptools import setup, find_packages

setup(
    name=\"fusarium_ma_pipeline\",
    version=\"1.0.0\",
    author=\"Manus AI\",
    description=\"Customized mutation accumulation pipeline for Fusarium graminearum\",
    long_description=open(\"README.md\").read(),
    long_description_content_type=\"text/markdown\",
    packages=find_packages(),
    scripts=[
        \n        \n        \"scripts/run_pipeline.py\",
        \"scripts/check_and_install_dependencies.sh\"
    ],
    install_requires=[
        \"numpy\",
        \"scipy\",
        \"pandas\",
        \"matplotlib\",
        \"seaborn\",
        \"pyyaml\",
        \"openpyxl\",
        \"click\"
    ],
    classifiers=[
        \"Programming Language :: Python :: 3\",
        \"License :: OSI Approved :: MIT License\",
        \"Operating System :: OS Independent\",
        \"Topic :: Scientific/Engineering :: Bio-Informatics\"
    ],
    python_requires=\">=3.8\",
)
