1.  Changing code to (not) throw an exception is considered
    **corrective**.

2.  Application code has precedence over test-code (i.e. the label ought
    to be determined by the application code first). If a e.g. a fix in
    app-code was made and unit tests were added (rule \#3, perfective),
    the label would still be corrective.

    A.  However, if the altered application code is considerably smaller
        in size and the rest of the changes are predominantly of a
        different change-type, then that label ought to be used.

3.  On very low confidence (do not forget to search for tickets using
    SHA1 and try keyword analysis), skip commit. This should probably
    also be applied to very large/tangled commits.

4.  If existing behavior is fixed to e.g. (not) cause an Exception/Error
    (usually w/o unit tests altered), then this is considered
    **corrective**. Also true for when the test(s) have not been altered
    but do pass now (i.e. non-passing test) or when tests have been
    altered to accomodate the corrected behavior.

    A.  If the corrective commit was made to *tighten* the code, then my
        rule \#7 need to be considered.

5.  Rule \#3 still applies if unit-tests were removed, if those removed
    parts can/should be considered now obsolete and are presumable
    superseded by the new test code.

6.  When **only** removing obsolete code resulting in required
    refactorings should be considered **perfective**.

    A.  Commits that merge functionality, without being actual
        merge-commits, should be considered to be **perfective**.

    B.  If existing code was **replaced** (subtle changes) to accomodate
        a new feature (e.g. adding a parameter in constructor and
        adapting the initialization of variables), then the change may
        be considered **adaptive**.

7.  *Tightening* of code, such as introducing generics or altering them
    or restricting the visibility of e.g. constructors, can be
    considered **perfective**, as it improves the code.

    A.  If a *seemingly* (i.e. not *outright*) corrective commit alters
        the code's behavior to *tighten* it (e.g. by now properly
        handling an NPE or by removing non-deterministic behavior), then
        this rule should be applied.

    B.  Improving of unit-tests should be considered **perfective**.

    C.  Modifications to code with concern to potential future bugs (and
        without a ticket) should also be considered **perfective**.

    D.  If the commit is *seemingly* tightening, but leaves the previous
        functionality in place (e.g. through overriding), then it may be
        considered **adapting**, as it add a new API.

8.  If new performance-tests were added, the commit is to be treated as
    if it was new unit tests (rule \#3).

9.  Always try to analyze the message using the manual commit
    classification helper, to see if there is a clear indication.

    A.  If all the scores are zero, attempt to replace words with
        synonyms that are actually matched, e.g. \"enable\" is a synonym
        for \"allow\".

10. If only renamings were made **and** the code cannot be said to be
    either corrective or adaptive, it should be **perfective**. In case
    of tangled commits, where my rule \#4 would normally apply, this
    rule has **precedence**.

11. **Loosening** or **exposing** of code in shape of removing/loosening
    of restrictions, such as exposing code by making it *public*, may be
    considered **adaptive**, as it \"adds\" the feature by exposing it,
    such that it appears to be a new feature.

12. Rule \#1 has low precedence, if the code changes were equal to or
    greater in size. In that case, the changes to the code need to be
    evaluated.

    A.  Updating *Javadoc* or *comments* of Rule \#1 shall refer only to
        code-accompanying comments, and **not** to other information to
        be found occasionally in code, such as license headers.

13. Refactorings, such as consolidation of functionalities, that do not
    correct bugs and do not introduce new functionality (except for
    functionality necessary to accomodate the refactoring), are to be
    considered **perfective**.

14. Reverting the changes made by a commit demands a proper
    investigation to its nature, especially if it reverts an already
    labelled commit (i.e. one cannot necessarily use the same label).
