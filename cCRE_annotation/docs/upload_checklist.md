# Prepare the review/upload package

From cCRE_annotation run:
~~~bash
python scripts/run_demo.py --outdir demo_output_review
python -m unittest discover -s tests -v
python scripts/verify_results.py
python scripts/make_submission.py
~~~

The archive contains source, versioned dependency instructions, demo data and
expected outputs, regression tests, documentation, MIT license, and available
global results. A SHA256SUMS.tsv inside the ZIP identifies every packaged file.
The archive is generated under dist/ and should not be added to normal git history.

For GitHub Desktop, review changes under cCRE_annotation, commit and push when
ready. Do not commit .venv, demo_output*, __pycache__, personal config/config.sh,
raw bigWigs or run intermediates. The provided .gitignore covers generated
directories. The optional developmental table may be omitted from a future
global-only release, but update the documentation/manifest consistently.

The existing result files are left in their supplied uncompressed format.
Do not silently replace them with a newly classified dataset after changing the
distance mode. Report a new analysis version and compare labels by ID.

Before submitting the Nature form, verify the remaining items in
nature_checklist.md: actual native/cluster versions and execution, manuscript
method alignment, original input access, reviewer access and independent
colleague installation. Do not claim those steps have been completed solely
because a form checkbox was previously selected.
