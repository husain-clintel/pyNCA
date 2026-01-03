"""Tests for parallel processing utilities."""

import numpy as np
import pytest

from pynca.utils.parallel import (
    get_n_workers,
    parallel_map,
    ParallelNCA,
)


class TestGetNWorkers:
    """Test get_n_workers function."""

    def test_all_cpus(self):
        """Test -1 returns all CPUs."""
        n = get_n_workers(-1)
        assert n >= 1

    def test_all_but_one(self):
        """Test -2 returns all but one CPU."""
        n = get_n_workers(-2)
        assert n >= 1

    def test_specific_number(self):
        """Test specific number of workers."""
        n = get_n_workers(2)
        assert n == 2 or n == get_n_workers(-1)  # May be capped at max

    def test_invalid_returns_one(self):
        """Test invalid input returns 1."""
        n = get_n_workers(0)
        assert n == 1


class TestParallelMap:
    """Test parallel_map function."""

    def test_basic(self):
        """Test basic parallel map."""
        def square(x):
            return x ** 2

        results = parallel_map(square, [1, 2, 3, 4, 5], n_jobs=2)
        assert results == [1, 4, 9, 16, 25]

    def test_single_item(self):
        """Test with single item (no parallelization)."""
        def double(x):
            return x * 2

        results = parallel_map(double, [5], n_jobs=2)
        assert results == [10]

    def test_empty_list(self):
        """Test with empty list."""
        def identity(x):
            return x

        results = parallel_map(identity, [], n_jobs=2)
        assert results == []

    def test_sequential_fallback(self):
        """Test fallback to sequential with n_jobs=1."""
        def cube(x):
            return x ** 3

        results = parallel_map(cube, [1, 2, 3], n_jobs=1)
        assert results == [1, 8, 27]


class TestParallelNCA:
    """Test ParallelNCA class."""

    def test_instantiation(self):
        """Test ParallelNCA instantiation."""
        pnca = ParallelNCA(n_jobs=2, backend="thread")
        assert pnca.n_jobs == 2
        assert pnca.backend == "thread"
        assert pnca.n_workers == 2

    def test_repr(self):
        """Test string representation."""
        pnca = ParallelNCA(n_jobs=2)
        repr_str = repr(pnca)
        assert "ParallelNCA" in repr_str
        assert "n_jobs=2" in repr_str

    def test_n_workers_property(self):
        """Test n_workers property."""
        pnca = ParallelNCA(n_jobs=-1)
        assert pnca.n_workers >= 1
