PYTHON ?= .venv/bin/python

.PHONY: setup doctor verify verify-rank2 integrity verify-all notes bundle
setup:
	bash scripts/setup.sh
doctor:
	$(PYTHON) scripts/doctor.py
verify:
	$(PYTHON) scripts/verify.py --suite quick
verify-rank2:
	$(PYTHON) scripts/verify.py --suite rank2
integrity:
	$(PYTHON) scripts/verify.py --suite integrity
verify-all:
	$(PYTHON) scripts/verify.py --suite all
notes:
	bash scripts/build_notes.sh
bundle:
	$(PYTHON) scripts/make_bundle.py
