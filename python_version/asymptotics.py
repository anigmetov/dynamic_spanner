#!/usr/bin/env python3

import numpy as np
import joblib as jl


def helper_prod(n, i):
    result = 1.0
    for j in range(i - 1):
        result = result * (1.0 - 1.0 / (n - j))
    result = i * result / (n - i)
    return result


def expectation(n):
    summands = jl.Parallel(n_jobs=-1)(jl.delayed(helper_prod)(n, i)
                                      for i in range(1, n))
    return sum(summands)
    # return sum([ helper_prod(n, i) for i in range(1, n)])


def optimistic_prob_matrix(n):
    p = np.zeros((n, n))
    for i in range(n):
        p[0][i] = 1.0 / n
    for r in range(1, n):
        for i in range(n):
            for j in range(i + 1, n):
                p[r][i] = p[r][i] + p[r - 1][j] * (1.0 / j)
    return p


def optimistic_expectation(n):
    p = optimistic_prob_matrix(n)
    p1 = p[:, 0]
    print(np.shape(p1))
    print(np.shape(np.array(range(n))))
    return np.dot(p1, np.array(range(n)))


def run_game(n):
    np.random.seed()
    if n == 1:
        return 0
    a = np.random.randint(low=1, high=n)
    rounds = 1
    while a != 1:
        a = np.random.randint(low=1, high=a)
        rounds = rounds + 1
    return rounds


def game_expectation(n, max_n_attempts=1000):
    game_results = jl.Parallel(n_jobs=-1)(jl.delayed(run_game)(n)
                                          for n_attempts in range(0,max_n_attempts))
    # game_results = [run_game(n) for attempt in range(0, max_n_attempts)]
    return sum(game_results) / len(game_results)


def markov_chain(s):
    n = s - 1
    Q = np.zeros((n, n))
    for row in range(n):
        for col in range(row+1, n):
            Q[row][col] = 1.0 / (n - row)
    N = np.linalg.inv(np.eye(n) - Q)
    print(N)
    print("********************************************************************************")
    # print(N[0])
    return (s, np.dot(N, np.ones((n, 1)))[0][0])


if __name__ == "__main__":
    # markov_chain(3)
    markov_chain(4)
    markov_chain(5)
    markov_chain(6)
    markov_chain(7)
    markov_chain(8)
    markov_chain(9)
    # results = jl.Parallel(n_jobs=-1)(jl.delayed(markov_chain)(2 ** s) for s in range(3,14))
    # for r in results:
    #     print(f"{r[0]} -> {r[1]}")
    # for n in range(25):
    #     # print(f"{2**n} -> {expectation(2 ** n)}")
    #     # print(f"{2**n} -> {optimistic_expectation(2 ** n)}")
    #     print(f"{n} -> {2**n} -> {game_expectation(2**n, min(2**n, 1000))}")
