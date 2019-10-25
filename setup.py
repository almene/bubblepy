from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='BubbleClustering',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Bubble clustering 3.7.3 python module with SNP profile import functionality",
    author="Amanda Saunders",
    author_email='saunders.mandy@hotmail.com',
    url='https://github.com/almene/BubbleClustering',
    packages=['bubblepy'],
    entry_points={
        'console_scripts': [
            'bubblepy=bubblepy.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='BubbleClustering',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ]
)
