#!/usr/bin/env python3
"""Reference grader for the two problems. Python 3 standard library only.

Usage:
    python3 -I grade.py p1 < model_response.txt
    python3 -I grade.py p2 < model_response.txt

Grades the terminal scientific answer only. Prose, rationale wording and formatting
outside the <EVAL_ANSWER> block are ignored by design.
Exit code 0 = PASS, 1 = FAIL.
"""
import json
import re
import sys

TAG = re.compile(r'<EVAL_ANSWER>(.*?)</EVAL_ANSWER>', re.S)


FIELDS = ['nonenhancing_fraction_treated', 'nonenhancing_fraction_placebo',
          'absolute_difference_pp', 'effect_on_margin_drug_exposure',
          'effect_on_margin_cytotoxic_sensitivity', 'recommended_endpoint_change',
          'answer']


def salvage(text):
    """Fallback parse when the <EVAL_ANSWER> block is absent.

    Added 2026-11-05 after LLM testing. Both models tested produced fully correct
    science while omitting the tags, and the brief instructs us to grade the terminal
    scientific answer rather than writing style; auto-failing a correct analysis on a
    missing wrapper would be the kind of formatting trick the brief rules out. So a
    missing block is recorded as a format deviation, not a scientific failure, provided
    every graded field can be recovered unambiguously from the prose.

    A field is only recovered when it appears exactly once, so an ambiguous response
    still fails rather than being guessed at.
    """
    out = {}
    for f in FIELDS:
        hits = re.findall(
            re.escape(f) + r'["\']?\s*[:=]\s*["\']?([A-Za-z_0-9.\-]+)', text)
        if len(set(hits)) == 1:
            v = hits[0]
            try:
                out[f] = float(v) if re.fullmatch(r'-?\d+(\.\d+)?', v) else v
            except ValueError:
                out[f] = v
    # bare "<N> percentage points" for the difference, and a final letter answer
    if 'absolute_difference_pp' not in out:
        m = re.findall(r'(\d+\.?\d*)\s*percentage points', text)
        if m and len(set(m)) >= 1:
            out['absolute_difference_pp'] = float(m[-1])
    if 'answer' not in out:
        m = re.findall(r'\b([A-D])\b(?=\s+is the best-supported|\s*is best supported)', text)
        if len(set(m)) == 1:
            out['answer'] = m[0]
    return out


def extract(text):
    """Prefer the structured block; fall back to prose recovery (see salvage)."""
    blocks = TAG.findall(text)
    if len(blocks) > 1:
        return None, None, 'found %d <EVAL_ANSWER> blocks, expected exactly 1' % len(blocks)
    if len(blocks) == 1:
        try:
            return json.loads(blocks[0].strip()), False, None
        except Exception as e:
            return None, None, 'EVAL_ANSWER block is not valid JSON: %s' % e
    d = salvage(text)
    if not d:
        return None, None, 'no <EVAL_ANSWER> block and no graded field recoverable from prose'
    return d, True, None


def num(d, k):
    v = d.get(k)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    return float(v)


def grade_p1(d):
    """Problem 1. Fractions may be given as proportions or as percents: both are
    unambiguous here, so both are accepted (documented in the grading rule)."""
    fails = []

    def frac(key, target):
        v = num(d, key)
        if v is None:
            fails.append('%s missing or not numeric' % key)
            return
        ok = abs(v - target) <= 0.01 or abs(v - target * 100) <= 1.0
        if not ok:
            fails.append('%s = %s, expected %.2f (or %.1f as percent)'
                         % (key, v, target, target * 100))

    frac('nonenhancing_fraction_treated', 0.4417)
    frac('nonenhancing_fraction_placebo', 0.1480)

    v = num(d, 'absolute_difference_pp')
    if v is None:
        fails.append('absolute_difference_pp missing or not numeric')
    elif not (abs(v - 29.4) <= 1.0 or abs(v - 0.294) <= 0.01):
        fails.append('absolute_difference_pp = %s, expected 29.4 pp' % v)

    a = str(d.get('answer', '')).strip().upper()
    if a != 'B':
        fails.append('answer = %r, expected B' % a)
    return fails


def grade_p2(d):
    """Problem 2. Categorical, exact match, case-insensitive.

    `unchanged` on exposure is accepted only when the rationale justifies it by the
    margin barrier already being intact. That judgment is not mechanisable, so the
    script flags it for a human rather than silently passing or failing it."""
    fails = []
    exp = str(d.get('effect_on_margin_drug_exposure', '')).strip().lower()
    if exp == 'unchanged':
        print('  [REVIEW] exposure = "unchanged": passes ONLY if the rationale argues '
              'the margin barrier was already intact. Check the rationale by hand.')
    elif exp != 'decreases':
        fails.append('effect_on_margin_drug_exposure = %r, expected decreases' % exp)

    sens = str(d.get('effect_on_margin_cytotoxic_sensitivity', '')).strip().lower()
    if sens != 'decreases':
        fails.append('effect_on_margin_cytotoxic_sensitivity = %r, expected decreases' % sens)

    ep = str(d.get('recommended_endpoint_change', '')).strip().lower()
    if ep != 'adopt_nonenhancing_inclusive_criteria':
        fails.append('recommended_endpoint_change = %r, expected '
                     'adopt_nonenhancing_inclusive_criteria' % ep)

    a = str(d.get('answer', '')).strip().upper()
    if a != 'B':
        fails.append('answer = %r, expected B' % a)
    return fails


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ('p1', 'p2'):
        sys.exit(__doc__)
    d, salvaged, err = extract(sys.stdin.read())
    if err:
        print('FAIL — %s' % err)
        sys.exit(1)
    missing = [f for f in (['nonenhancing_fraction_treated',
                            'nonenhancing_fraction_placebo',
                            'absolute_difference_pp', 'answer']
                           if sys.argv[1] == 'p1' else
                           ['effect_on_margin_drug_exposure',
                            'effect_on_margin_cytotoxic_sensitivity',
                            'recommended_endpoint_change', 'answer'])
               if f not in d]
    if missing:
        print('FAIL — graded field(s) not recoverable: %s' % ', '.join(missing))
        sys.exit(1)
    fails = grade_p1(d) if sys.argv[1] == 'p1' else grade_p2(d)
    if fails:
        print('FAIL')
        for f in fails:
            print('  - %s' % f)
        sys.exit(1)
    if salvaged:
        print('PASS (science correct; FORMAT DEVIATION — no <EVAL_ANSWER> block, '
              'fields recovered from prose)')
    else:
        print('PASS')
    sys.exit(0)


if __name__ == '__main__':
    main()
