from setuptools import setup

setup(
    name = 'gcmnester',
    version = '0.1',
    description = "A tool for building nested models in MITgcm.", 
    url = 'http://github.com/glwagner/gcmnester',
    author = 'Gregory L. Wagner',
    author_email = 'wagner.greg@gmail.com',
    license = 'MIT',
    packages = ['gcmnester'],
    install_requires = [
        'numpy', 
        'matplotlib', 
        'globotopo',
    ],
    zip_safe = False,
)
