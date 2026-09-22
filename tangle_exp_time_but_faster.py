#AI-gen
from __future__ import annotations

import argparse
import itertools
import json
import time
from collections import Counter, defaultdict
from functools import lru_cache
from typing import Iterable


Tree = tuple[int, int]
Perm = tuple[int, ...]  # zero-based image notation: p[left_leaf] = right_leaf
Tanglegram = tuple[Tree, Tree, Perm]


@lru_cache(maxsize=None)
def W(n: int) -> int:
    """Wedderburn-Etherington number for n leaves."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    if n % 2:
        return sum(W(i) * W(n - i) for i in range(1, (n + 1) // 2))
    half = n // 2
    return (
        sum(W(i) * W(n - i) for i in range(1, half))
        + W(half) * (W(half) + 1) // 2
    )


def root_split_finder(pair: Tree) -> tuple[Tree, Tree]:
    """Return the root split of the s-th unordered binary tree on n leaves."""
    n, t = pair
    if t >= W(n):
        raise ValueError(f"there is no tree for t={t} and n={n}")
    if n == 1:
        return (1, 0), (0, 0)

    if n % 2:
        for i in range(1, (n + 1) // 2):
            block = W(n - i) * W(i)
            if t < block:
                s_big, s_small = divmod(t, W(i))
                return (n - i, s_big), (i, s_small)
            t -= block
    else:
        half = n // 2
        for i in range(1, half):
            block = W(n - i) * W(i)
            if t < block:
                s_big, s_small = divmod(t, W(i))
                return (n - i, s_big), (i, s_small)
            t -= block

        y = W(half)
        for small_index in range(y):
            width = y - small_index
            if t < width:
                return (half, small_index + t), (half, small_index)
            t -= width

    raise AssertionError("unreachable tree index")


def binary_tree(tree_0: Tree, tree_1: Tree) -> Tree:
    """Return the unordered binary tree tree_0 + tree_1 in WE order."""
    (n_0, s_0), (n_1, s_1) = sorted((tree_0, tree_1))
    if n_0 == 0:
        return n_1, s_1
    if n_1 == 0:
        return n_0, s_0

    n = n_0 + n_1
    offset = sum(W(n - i) * W(i) for i in range(1, n_0))
    if n_0 == n_1:
        y = W(n_0)
        return n, offset + sum(y - i for i in range(s_0)) + s_1 - s_0
    return n, offset + W(n_0) * s_1 + s_0


@lru_cache(maxsize=None)
def tree_automorphisms(tree: Tree) -> tuple[Perm, ...]:
    """
    Automorphisms of a rooted unordered binary tree, acting on its leaves.

    A permutation a is represented as a tuple with a[old_leaf] = new_leaf.
    """
    n, _ = tree
    if n <= 1:
        return (tuple(range(n)),)

    left, right = root_split_finder(tree)
    n_left = left[0]
    left_autos = tree_automorphisms(left)
    right_autos = tree_automorphisms(right)

    autos: list[Perm] = []
    for a_left in left_autos:
        for a_right in right_autos:
            image = [0] * n
            for i, v in enumerate(a_left):
                image[i] = v
            for i, v in enumerate(a_right):
                image[n_left + i] = n_left + v
            autos.append(tuple(image))

            if left == right:
                swapped = [0] * n
                for i, v in enumerate(a_left):
                    swapped[i] = n_left + v
                for i, v in enumerate(a_right):
                    swapped[n_left + i] = v
                autos.append(tuple(swapped))

    return tuple(sorted(set(autos)))


@lru_cache(maxsize=None)
def canonical_tuple_under_tree_aut(tree: Tree, tup: tuple[int, ...]) -> tuple[int, ...]:
    """
    Smallest tuple obtainable by applying automorphisms of tree to positions.

    This one-sided minimization is correct and fast. The unsound step was trying
    to do both tanglegram sides independently in sequence.
    """
    n, _ = tree
    if n <= 1:
        return tup

    left, right = root_split_finder(tree)
    n_left = left[0]
    tup_left = canonical_tuple_under_tree_aut(left, tup[:n_left])
    tup_right = canonical_tuple_under_tree_aut(right, tup[n_left:])

    unswapped = tup_left + tup_right
    if left == right:
        return min(unswapped, tup_right + tup_left)
    return unswapped


def inverse_perm(p: Perm) -> Perm:
    inv = [0] * len(p)
    for i, v in enumerate(p):
        inv[v] = i
    return tuple(inv)


@lru_cache(maxsize=None)
def canonical_values_under_tree_aut(tree: Tree, p: Perm) -> Perm:
    """Smallest tuple obtainable by applying automorphisms of tree to values."""
    return inverse_perm(canonical_tuple_under_tree_aut(tree, inverse_perm(p)))


@lru_cache(maxsize=None)
def canonical_perm(left_tree: Tree, right_tree: Tree, p: Perm) -> Perm:
    """
    Smallest permutation in the tanglegram isomorphism orbit.

    Exact hybrid:
      - enumerate automorphisms on the smaller side;
      - recursively minimize the other side.

    The fully staged right_canon(left_canon(p)) version is not exact.
    """
    left_autos = tree_automorphisms(left_tree)
    right_autos = tree_automorphisms(right_tree)
    best: Perm | None = None

    if len(left_autos) <= len(right_autos):
        for a_left in left_autos:
            q = [0] * len(p)
            for old_left, old_right in enumerate(p):
                q[a_left[old_left]] = old_right
            candidate = canonical_values_under_tree_aut(right_tree, tuple(q))
            if best is None or candidate < best:
                best = candidate
    else:
        for a_right in right_autos:
            q = tuple(a_right[value] for value in p)
            candidate = canonical_tuple_under_tree_aut(left_tree, q)
            if best is None or candidate < best:
                best = candidate

    assert best is not None
    return best


def canonical_tanglegram(tan: Tanglegram) -> Tanglegram:
    left_tree, right_tree, p = tan
    return left_tree, right_tree, canonical_perm(left_tree, right_tree, tuple(p))


def is_isomorphic(tan_1: Tanglegram, tan_2: Tanglegram) -> bool:
    if tan_1[0] != tan_2[0] or tan_1[1] != tan_2[1]:
        return False
    return canonical_tanglegram(tan_1) == canonical_tanglegram(tan_2)


def mid(button: tuple[int, int]) -> int:
    a, b = button
    if a == b:
        return a - 1
    return (a + b - 1) // 2


def press_buttons(p: Perm, dct: dict[int, list[tuple[int, int]]]) -> Perm:
    """Port of the Sage Permutation button logic, using zero-based tuples."""

    def press_left_buttons(q: Perm, buttons: list[tuple[int, int]]) -> Perm:
        for button in buttons:
            if not button:
                continue
            a, b = button
            m = mid(button)
            q = q[: a - 1] + q[m:b] + q[a - 1 : m] + q[b:]
        return q

    def press_right_buttons(q: Perm, buttons: list[tuple[int, int]]) -> Perm:
        for button in buttons:
            if not button:
                continue
            q_inv = inverse_perm(q)
            a, b = button
            m = mid(button)
            q_inv = q_inv[: a - 1] + q_inv[m:b] + q_inv[a - 1 : m] + q_inv[b:]
            q = inverse_perm(q_inv)
        return q

    return press_right_buttons(press_left_buttons(p, dct[0]), dct[1])


def canonicalize_after_deleting(values: Iterable[int], removed: int) -> Perm:
    """Relabel 0..n-1 with one value removed to 0..n-2, preserving order."""
    return tuple(x if x < removed else x - 1 for x in values)


def subtree_button(global_offset: int, n_leaves: int, removed_global: int) -> tuple[int, int]:
    """
    Return the 1-based first/last button for a subtree in the card permutation.

    global_offset is the original 0-based first leaf of this subtree.
    n_leaves includes the deleted leaf if this subtree contains it.
    removed_global is the original 0-based deleted leaf.
    """
    g, n, r = global_offset, n_leaves, removed_global
    if n <= 1:
        return (0, 0)
    if g + n - 1 < r:
        return (g + 1, g + n)
    if g > r:
        return (g, g + n - 1)
    return (g + 1, g + n - 1)


def tan_multideck_counter(tan: Tanglegram) -> Counter[Tanglegram]:
    """Return the multideck as a Counter keyed by canonical card tanglegrams."""
    cbus: defaultdict[int, dict[int, list[tuple[int, int]]]] = defaultdict(
        lambda: {0: [], 1: []}
    )

    def card_finder_with_left_buttons(
        pair: Tree, tup: tuple[int, ...], j: int, global_offset: int = 0
    ) -> Tree:
        n, _ = pair
        k = len(tup)
        removed_global = j - 1
        if k == 0:
            return 0, 0
        if k == n:
            return pair

        left_tree, right_tree = root_split_finder(pair)
        split = left_tree[0]
        left_sublist = tuple(i for i in tup if i < split)
        right_sublist = tuple(i for i in tup if i >= split)
        right_adjusted = tuple(i - split for i in right_sublist)

        left_card = card_finder_with_left_buttons(
            left_tree, left_sublist, j, global_offset
        )
        right_card = card_finder_with_left_buttons(
            right_tree, right_adjusted, j, global_offset + split
        )
        card = binary_tree(left_card, right_card)

        if left_card < right_card and 1 < split and k > 2:
            cbus[j][0].append(subtree_button(global_offset, n, removed_global))
        return card

    def card_finder_with_right_buttons(
        pair: Tree,
        tup: tuple[int, ...],
        deleted_right_label: int,
        p_inv: Perm,
        global_offset: int = 0,
    ) -> Tree:
        n, _ = pair
        k = len(tup)
        removed_global = deleted_right_label - 1
        if k == 0:
            return 0, 0
        if k == n:
            return pair

        left_tree, right_tree = root_split_finder(pair)
        split = left_tree[0]
        left_sublist = tuple(i for i in tup if i < split)
        right_sublist = tuple(i for i in tup if i >= split)
        right_adjusted = tuple(i - split for i in right_sublist)

        left_card = card_finder_with_right_buttons(
            left_tree, left_sublist, deleted_right_label, p_inv, global_offset
        )
        right_card = card_finder_with_right_buttons(
            right_tree,
            right_adjusted,
            deleted_right_label,
            p_inv,
            global_offset + split,
        )
        card = binary_tree(left_card, right_card)

        if left_card < right_card and 1 < split and k > 2:
            cbus[p_inv[deleted_right_label - 1] + 1][1].append(
                subtree_button(global_offset, n, removed_global)
            )
        return card

    left_tree, right_tree, p = tan
    n = left_tree[0]
    p_inv = inverse_perm(p)
    cards: Counter[Tanglegram] = Counter()

    for j in range(1, n + 1):
        deleted_right = p[j - 1] + 1
        subset_left = tuple(range(j - 1)) + tuple(range(j, n))
        subset_right_unsorted = tuple(p[k] for k in range(n) if k != j - 1)
        sorted_subset_right = tuple(range(deleted_right - 1)) + tuple(
            range(deleted_right, n)
        )

        card_left_tree = card_finder_with_left_buttons(left_tree, subset_left, j, 0)
        card_right_tree = card_finder_with_right_buttons(
            right_tree, sorted_subset_right, deleted_right, p_inv, 0
        )
        card_perm = canonicalize_after_deleting(
            subset_right_unsorted, deleted_right - 1
        )
        raw_card = (
            card_left_tree,
            card_right_tree,
            press_buttons(card_perm, cbus[j]),
        )
        cards[canonical_tanglegram(raw_card)] += 1

    return cards


def tan_deck_signature(tan: Tanglegram) -> frozenset[Tanglegram]:
    return frozenset(tan_multideck_counter(tan))


def generate_list_of_representatives(n: int) -> list[Tanglegram]:
    """Generate one canonical representative for each tanglegram class."""
    representatives: dict[Tanglegram, Tanglegram] = {}
    trees = [(n, s) for s in range(W(n))]
    perms = list(itertools.permutations(range(n)))

    for left_tree in trees:
        for right_tree in trees:
            for p in perms:
                canon = canonical_tanglegram((left_tree, right_tree, p))
                representatives.setdefault(canon, canon)

    return sorted(representatives)


def tanglegram_for_display(tan: Tanglegram) -> tuple[Tree, Tree, list[int]]:
    left_tree, right_tree, p = tan
    return left_tree, right_tree, [x + 1 for x in p]


def tanglegram_as_sageish(tan: Tanglegram) -> str:
    left_tree, right_tree, p = tanglegram_for_display(tan)
    return f"({left_tree}, {right_tree}, Permutation({p}))"


def counter_as_jsonable(counter: Counter[Tanglegram]) -> list[dict[str, object]]:
    return [
        {"card": tanglegram_for_display(card), "count": count}
        for card, count in sorted(counter.items())
    ]


def find_counterexamples(
    n: int, output: str | None = None, limit: int | None = None
) -> list[dict[str, object]]:
    started = time.perf_counter()
    reps = generate_list_of_representatives(n)
    print(f"Found {len(reps)} unique tanglegram classes for n={n}.")

    deck_groups: defaultdict[frozenset[Tanglegram], list[int]] = defaultdict(list)
    multidecks: list[Counter[Tanglegram]] = []
    for idx, tan in enumerate(reps, start=1):
        multideck = tan_multideck_counter(tan)
        multidecks.append(multideck)
        deck_groups[frozenset(multideck)].append(idx - 1)
        if idx % 1000 == 0:
            elapsed = time.perf_counter() - started
            print(f"  computed decks for {idx}/{len(reps)} reps in {elapsed:.1f}s")

    candidate_groups = [indices for indices in deck_groups.values() if len(indices) > 1]
    candidate_pair_count = sum(len(g) * (len(g) - 1) // 2 for g in candidate_groups)
    print(
        f"Finished decks: {len(deck_groups)} deck signatures, "
        f"{candidate_pair_count} equal-deck candidate pairs."
    )

    results: list[dict[str, object]] = []
    out_handle = open(output, "w", encoding="utf-8") if output else None
    try:
        for group in candidate_groups:
            for pos, i in enumerate(group):
                for j in group[pos + 1 :]:
                    rep_1 = reps[i]
                    rep_2 = reps[j]
                    multidecks_equal = multidecks[i] == multidecks[j]
                    inverse_switch = inverse_perm(rep_1[2]) == rep_2[2]
                    record = {
                        "rep_1": tanglegram_for_display(rep_1),
                        "rep_2": tanglegram_for_display(rep_2),
                        "multidecks_equal": multidecks_equal,
                        "inverse_switch": inverse_switch,
                    }
                    results.append(record)

                    if out_handle:
                        out_handle.write(json.dumps(record) + "\n")

                    if limit is not None and len(results) >= limit:
                        print(f"Stopped after --limit={limit} counterexample pairs.")
                        return results
    finally:
        if out_handle:
            out_handle.close()

    elapsed = time.perf_counter() - started
    print(f"Found {len(results)} counterexample pairs in {elapsed:.1f}s.")
    return results


def self_test() -> None:
    expected_classes = {
        1: 1,
        2: 1,
        3: 2,
        4: 13,
        5: 114,
        6: 1509,
        7: 25595,
    }
    for n, expected in expected_classes.items():
        got = len(generate_list_of_representatives(n))
        status = "ok" if got == expected else "FAIL"
        print(f"classes n={n}: got {got}, expected {expected} [{status}]")
        if got != expected:
            raise AssertionError((n, got, expected))


def _from_display(display: object) -> Tanglegram:
    left, right, perm = display  # type: ignore[misc]
    return tuple(left), tuple(right), tuple(x - 1 for x in perm)  # type: ignore[arg-type]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fast tanglegram deck search with exact hybrid canonicalization."
    )
    parser.add_argument("-n", type=int, default=7)
    parser.add_argument("--output", help="Optional JSONL file for counterexample pairs.")
    parser.add_argument(
        "--limit", type=int, help="Stop after this many counterexample pairs."
    )
    parser.add_argument(
        "--self-test", action="store_true", help="Check known tanglegram class counts."
    )
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    results = find_counterexamples(args.n, args.output, args.limit)
    if results:
        print("First pair:")
        print("  ", tanglegram_as_sageish(_from_display(results[0]["rep_1"])))
        print("  ", tanglegram_as_sageish(_from_display(results[0]["rep_2"])))


if __name__ == "__main__":
    main()
