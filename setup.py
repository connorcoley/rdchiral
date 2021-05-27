import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rdchiral",
    version="1.1.0",
    author="Connor Coley",
    author_email="ccoley@mit.edu",
    description="Wrapper for RDKit's RunReactants to improve stereochemistry handling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/connorcoley/rdchiral/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)
