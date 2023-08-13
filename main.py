from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dwave.inspector
from collections import defaultdict


qubo_name = "H"

strings = ["aaa", "aaa", "ddd"]
# strings = ["aaa", "aaa", "ddd", "ddd", "ddd"]
# strings = ["aaa", "aaa", "ded", "ded", "ded", "ddd"]
# strings = [
#     "abcdef",
#     "ghijkl",
#     "abcghi",
#     "xyzjkl",
#     "abcmno"
# ]

A = 3
B = 1


def f(c1, c2):
    return 1 if c1 != c2 else 0


def main():
    Q = defaultdict(int)

    n = len(strings)
    m = len(strings[0])

    for i in range(m):
        offset = i * n
        for j in range(n):
            if qubo_name == "H":
                Q[(offset + j, offset + j)] = B * sum([f(string[i], strings[j][i]) for string in strings]) - A
            elif qubo_name == "H'":
                Q[(offset + j, offset + j)] = B * sum([((ord(strings[j][i]) - ord(string[i])) ** 2) / (((ord(strings[j][i]) - ord(string[i])) ** 2) + 1) for string in strings]) - A

            for k in range(j + 1, n):
                Q[(offset + j, offset + k)] = A

    # ------- Run our QUBO on the QPU -------
    # Set up QPU parameters
    chainstrength = 8
    numruns = 100

    sampler = EmbeddingComposite(DWaveSampler())
    response = sampler.sample_qubo(Q,
                                   chain_strength=chainstrength,
                                   num_reads=numruns,
                                   label='Example - Closest String')

    outputs = {}
    total_occurrences = 0

    print('-' * 130)
    print('{:>70s}{:>30s}{:^30s}'.format('Raw Output', 'String', 'Energy'))
    print('-' * 130)
    for sample, E, num_occurrences in response.data(fields=['sample', 'energy', 'num_occurrences']):
        flat_one_hot_clst_str = [v for k, v in sample.items()]
        assert len(flat_one_hot_clst_str) == m * n

        output = ""
        for i in range(m):
            offset = i * n

            found_char = False
            for j in range(n):
                found_char = flat_one_hot_clst_str[offset + j] == 1
                if found_char:
                    output += strings[j][i]
                    break

            if not found_char:
                output += '_'

        if output not in outputs:
            outputs[output] = {
                "raw_outputs": [flat_one_hot_clst_str],
                "cumulative_energy": E,
                "num_occurrences": num_occurrences,
                "count": 1
            }
        else:
            outputs[output]["raw_outputs"].append(flat_one_hot_clst_str)
            outputs[output]["cumulative_energy"] += E
            outputs[output]["num_occurrences"] += num_occurrences
            outputs[output]["count"] += 1

        total_occurrences += num_occurrences

        print('{:>70s}{:>30s}{:^30s}'.format(str(flat_one_hot_clst_str), output, str(E)))

    print()

    print('-' * 130)
    print('{:>30s}{:>30s}{:^30s}'.format('Output', 'Average Energy', 'Occurrence Ratio'))
    print('-' * 130)
    for output, output_data in outputs.items():
        print('{:>30s}{:>30f}{:^30f}'.format(output, output_data["cumulative_energy"] / output_data["count"], output_data["num_occurrences"] / total_occurrences))

    dwave.inspector.show(response)


if __name__ == "__main__":
    main()
