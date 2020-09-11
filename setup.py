from setuptools import setup, find_packages
import os

module_dir = os.path.dirname(os.path.abspath(__file__))
reqs_raw = open(os.path.join(module_dir, "requirements.txt")).read()
# reqs_list = [r.replace("==", ">=") for r in reqs_raw.split("\n")]
reqs_list = [r for r in reqs_raw.split("\n")]


if __name__ == "__main__":

    with open(os.path.join(module_dir, "VERSION"), "r") as f:
        version = f.read()

    setup(
        name='amlt',
        version=version,
        description='automated generation of MLIP training sets',
        long_description="Automated generation of machine learned interatomic potential training sets"
                         "https://github.com/MTD-group/amlt",
        url='https://github.com/MTD-group/amlt',
        author=['Michael Waters', 'Nicholas Wagner'],
        author_email='michael.j.waters@northwestern.edu',
        license='modified BSD',
        packages=find_packages(where=".", exclude=("benchdev", "benchdev.*")),
        package_data={},
        zip_safe=False,
        install_requires=reqs_list,
        extras_require={},
        classifiers=['Programming Language :: Python :: 3.6',
                     'Development Status :: 3 - Alpha',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Intended Audience :: Information Technology',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        test_suite='amlt',
        tests_require='tests'
    )
