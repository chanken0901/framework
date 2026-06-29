# Project-level Makefile

.PHONY: help case-list clean

help:
	@echo "Available targets:"
	@echo "  make case-list   Show registered cases"
	@echo "  make clean       Remove common temporary files"

case-list:
	@if [ -f cases/case_index.csv ]; then 		cat cases/case_index.csv; 	else 		echo "cases/case_index.csv not found"; 	fi

clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
