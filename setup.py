import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='genomeencode',
    version='1.0.0',
    description='Optimized DNA Sequence Compression and Reconstruction using Burros-Wheeler and Huffman coding algorithms',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Bhavanapoli08/genomeencode',
    author='Bhavana Poli',
    author_email = 'bhavanapoli61@gmail.com',
    download_url = 'https://github.com/Bhavanapoli08/genomeencode/archive/v1.0.0.tar.gz',
    license='MIT',
    packages=setuptools.find_packages(include=['genomeencode', 'genomeencode.*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #install_requires=REQUIREMENTS,
    extras_require={
        'dev': [
            'pytest >= 6.0.0',
            'pytest-cov >= 2.10.0',
            'coveralls >= 2.1.2',
            'flake8 >= 3.8.0',
            'mock >= 4.0.0',
        ]
    },
    entry_points={
        'gui_scripts': [
            'genomeencode=genomeencode.interface:Interface.main'
        ]
    },
    python_requires='>=3.6',
)
