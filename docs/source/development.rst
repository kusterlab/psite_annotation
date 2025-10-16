Development
===========

To add a new annotator to the package, follow these steps:

#. Make a fork of the repository.

   - Go to the project's GitHub page (https://github.com/kusterlab/psite_annotation).
   - In the top-right corner of the page, click Fork.
   - Choose your personal GitHub account (or organization) as the destination.

#. Make a clone of your forked repository:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/psite_annotation.git

#. Create a new file in `psite_annotation/annotators`, the easiest is to make a copy of an existing annotator such as `domain.py`.
#. Implement the `__init__()`, `load_annotations()` and `annotate()` functions.
#. Create unit tests in `tests/unit_tests/annotators`, strive for 100% code coverage and covering edge cases.
#. Add a function that uses your annotator to `psite_annotation/functional_annotation.py` and add that function's name to the `__all__` list in the top of that file.
#. Create a system test for that function in `tests/system_tests/test_functional_annotation.py`.
#. Create a pull request:

   - Go to your fork on GitHub (e.g. https://github.com/YOUR_USERNAME/psite_annotation).
   - You should see a “Compare & pull request” button — click it. Alter=natively, go to the Pull Requests tab and click “New pull request”.