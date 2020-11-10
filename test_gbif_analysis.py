#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
description
"""

import gbif_analysis
import pandas as pd
import numpy as np


def test_cooccurrence_matrix():
    """
    Tests build_cooccurrence_matrix() with a simple DataFrame containing 0 and 1.
    """
    # Given
    df = pd.DataFrame({'Basket': ['A', 'B', 'C', 'D'],
                       'Apples': [1, 0, 1, 1],
                       'Bananas': [1, 1, 1, 0],
                       'Oranges': [1, 1, 1, 0]}).set_index('Basket')

    # When
    expected_cooc_m = pd.DataFrame([[3, 2, 2],
                                    [2, 3, 3],
                                    [2, 3, 3]],
                                   columns=['Apples', 'Bananas', 'Oranges'], index=['Apples', 'Bananas', 'Oranges'])
    expected_cooc_m_f = pd.DataFrame([[0, 2, 2],
                                      [2, 0, 3],
                                      [2, 3, 0]],
                                     columns=['Apples', 'Bananas', 'Oranges'], index=['Apples', 'Bananas', 'Oranges'])
    cooc_m, cooc_m_f = gbif_analysis.build_cooccurrence_matrix(df)

    # Then
    pd.testing.assert_frame_equal(cooc_m, expected_cooc_m, check_dtype=False)
    pd.testing.assert_frame_equal(cooc_m_f, expected_cooc_m_f, check_dtype=False)


def test_cooccurrence_matrix_with_nulls():
    """
    Tests build_cooccurrence_matrix() with a DataFrame containing integers and NaN-values.
    """
    # Given
    df = pd.DataFrame({'Basket': ['A', 'B', 'C', 'D'],
                       'Apples': [20, np.nan, 8, 13],
                       'Bananas': [3, 5, 1, 0],
                       'Oranges': [3, 7, 4, np.nan]}).set_index('Basket')

    # When
    expected_cooc_m = pd.DataFrame([[3, 2, 2],
                                    [2, 3, 3],
                                    [2, 3, 3]],
                                   columns=['Apples', 'Bananas', 'Oranges'], index=['Apples', 'Bananas', 'Oranges'])
    expected_cooc_m_f = pd.DataFrame([[0, 2, 2],
                                      [2, 0, 3],
                                      [2, 3, 0]],
                                     columns=['Apples', 'Bananas', 'Oranges'], index=['Apples', 'Bananas', 'Oranges'])
    cooc_m, cooc_m_f = gbif_analysis.build_cooccurrence_matrix(df)

    # Then
    pd.testing.assert_frame_equal(cooc_m, expected_cooc_m, check_dtype=False)
    pd.testing.assert_frame_equal(cooc_m_f, expected_cooc_m_f, check_dtype=False)
