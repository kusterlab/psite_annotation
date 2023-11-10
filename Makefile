package_build:
	poetry build

publish_test:
	poetry publish -r testpypi

publish:
	poetry publish

test: unit_test system_test

unit_test:
	python3 -m pytest --cov=psite_annotation --cov-report html --cov-report term -s tests/unit_tests

system_test:
	python3 -m pytest --cov=psite_annotation --cov-report html --cov-report term -s tests/system_tests
