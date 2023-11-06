test: unit_test system_test

unit_test:
	python3 -m pytest --cov=psite_annotation --cov-report html --cov-report term tests/unit_tests

system_test:
	python3 -m pytest --cov=psite_annotation --cov-report html --cov-report term tests/system_tests
