from setuptools import setup, find_packages

setup(
    name="lipidpy",
    version="0.1.0",
    description="Command line tool for untargeted lipid identification",
    author="biryb",
    author_email="birgittaryback@gmail.com",
    license="MIT",
    long_description_content_type="text/markdown",  # Assumes README is markdown
    long_description=open("README.md").read(),
    python_requires=">=3.11",
    install_requires=[
        "polars>=1.21.0,<2.0.0",
        "pyteomics>=4.7.5,<5.0.0",
        "fastexcel>=0.12.1,<0.13.0",
        "numpy>=2.2.2,<3.0.0"
    ],
    packages=find_packages(where="LipIDpy"),
    entry_points={
        "console_scripts": [
            "lipidpy=LipIDpy.lipidpy:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
