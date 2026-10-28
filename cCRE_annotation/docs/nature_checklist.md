# Mapping to the supplied Nature software checklist

Basis: the user-supplied nr-software-policy-20261009.pdf, whose internal form is
titled "Code and Software Submission Checklist" and carries "June 2017".
The filename date is not assumed to be the policy revision date. This document
maps that attachment, not a claim of editorial approval or all Nature policies.

| Attachment item | Evidence in this module | Status |
| --- | --- | --- |
| Source code | scripts/ | Present |
| Small real/simulated demonstration dataset | demo/input and demo/expected | Present; core demo verified |
| README software/OS dependencies and versions | README sections 1-2, requirements.txt, environment.yml | Core tested; Linux recipe explicitly unverified |
| Tested software versions | docs/validation.md | Core environment recorded |
| Non-standard hardware | README section 1 | None required |
| Installation instructions | README section 2 | Provided |
| Typical installation time | README and validation record | Estimate plus measured clean-install result |
| Demo instructions | README section 3 | Executed |
| Expected output | demo/expected and demo/README.md | Automatically checked |
| Expected demo runtime | README and demo report | Measured on this host; desktop estimate included |
| Instructions for own data | README section 4 and config templates | Provided; native full run pending |
| Optional quantitative reproduction instructions | verify_results.py and method_summary.md | Archived classifications/counts checked; raw rerun needs original inputs |
| License | LICENSE copied from parent repository | MIT, existing copyright retained |
| Repository link | README section 6 | Provided; reviewer access must be verified |
| Detailed functionality/pseudocode in manuscript | docs/method_summary.md | Ready for Methods alignment; manuscript itself not edited |
| Single ZIP or accessible review link | make_submission.py | ZIP generation supported |
| Unfamiliar colleague installation feedback | Not replaced by automated tests | Author action still required |

The PDF contains checked form fields; those are author declarations, not evidence
that the previous software package met every requirement. The filled form
identifies the manuscript Methods subsection "Identification of candidate
cis-regulatory elements in the macaque genome". Align that subsection with the
documented normalization and explicit TSS-distance behavior before submission.

Still to verify outside this Windows core validation:
- Native Linux bigWig extraction, bedtools invocation and chosen SLURM version.
- Optional R plot and versions on the production system.
- Original data references/accessions and exact reference/input hashes.
- Manuscript consistency if corrected distance mode is used.
- Installation by a colleague unfamiliar with the code.
- Public/reviewer access to the submitted version; this work does not push or release.

An immutable release tag/archive is useful for identifying reviewed code.
A DOI is not presented here as a requirement of the supplied checklist.
