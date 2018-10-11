#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>

const unsigned int GENOME_SIZE = 100;
const unsigned int POPULATION_SIZE = 20;
const unsigned int GENERATION_NUM = 5;
const unsigned int ELITE_POPULATION_SIZE = 4;
const double MUTATION_RATE = 0.00;

// 0 : Two-point Crossover, 1 : Uniform Crossover
const unsigned int TABLE = 1;

const double TABLE_INIT_COST = 1;

// 0 : Temporary, 1 : LifeGame (my function)
const unsigned int EVALUATION = 1;
// 0 : Temporary, 1 : ranking
const unsigned int SELECTION = 1;
// 0 : Temporary, 1 : Use crossover table
const unsigned int CROSSOVER = 0;
// 0 : Temporary
const unsigned int MUTATION = 0;

class Connection {
public:
	unsigned int node1, node2;
	double cost;
	Connection(unsigned int n1, unsigned int n2, double c) {
		node1 = n1;
		node2 = n2;
		cost = c;
		//std::cout << "(" << n1 << "," << n2 << "," << cost << ")" << std::endl;
	}
};

typedef std::vector<Connection> connections_t;
typedef std::vector<bool> genome_t;
typedef int fitness_t;

class Individual {
public:
	genome_t genome;
	fitness_t fitness;
	connections_t table, choice, cut;
	Individual() {}
	Individual(unsigned int genomeSize) {
		for (unsigned int i = 0; i < genomeSize; i++) {
			genome.insert(genome.begin(), rand() % 2);

			for (unsigned int j = i + 1; j < genomeSize; j++) {
				if (j == i + 1 || (i == 0 && j == genomeSize - 1)) {
					table.push_back(Connection(i, j, TABLE_INIT_COST));
				}
				else {
					if (TABLE == 0) {
						table.push_back(Connection(i, j, 0));
					}
					else if (TABLE == 1) {
						table.push_back(Connection(i, j, TABLE_INIT_COST));
					}
				}
			}
		}
	}
};

typedef std::vector<Individual> population_t;

class Population : public population_t {
public:
	Population() {}
	Population(unsigned int populationSize, unsigned int genomeSize) {
		for (unsigned int i = 0; i < populationSize; i++) {
			this->push_back(Individual(genomeSize));
		}
	}
	void show() {
		for (population_t::iterator it = population_t::begin(); it != population_t::end(); it++) {
			for (unsigned int j = 0; j < it->genome.size(); j++) {
				std::cout << it->genome[j];
			}
			std::cout << " " << it->fitness << std::endl;
		}
	}
};

typedef Population(*funcPtr_t)(Population);
typedef std::vector<funcPtr_t> funcPtrVec_t;

namespace Evaluation {
	Population temp(Population arg) {
		fitness_t count;
		for (unsigned int i = 0; i < arg.size(); i++) {
			count = 0;
			for (unsigned int j = 0; j < arg[i].genome.size(); j++) {
				if (arg[i].genome[j] == true) count++;
			}
			//std::cout << count << std::endl;
			arg[i].fitness = count;
		}
		return arg;
	}

	Population lifeGame(Population arg) {
		int size = ceil(sqrt(arg[0].genome.size()));
		unsigned int count;

		for (int i = 0; i < (int)arg.size(); i++) {
			arg[i].fitness = 0;
			for (int j = 0; j < (int)arg[i].genome.size(); j++) {
				count = 0;
				for (int k = -1; k <= 1; k++) {
					for (int l = -1; l <= 1; l++) {
						if (0 <= j + k + l*size && j + k + l*size < (int)arg[i].genome.size()) {
							if (arg[i].genome[j + k + l*size] == true) count++;
						}
					}
				}
				if (arg[i].genome[j] == true) {
					if (2 < count && count < 5) {
						arg[i].fitness += 1;
					}
					else {
						arg[i].fitness -= 1;
					}
				}
			}
		}
		return arg;
	}

	const funcPtrVec_t call = { temp, lifeGame };
}

const funcPtr_t evaluation = Evaluation::call[EVALUATION];

namespace Selection {
	Population temp(Population arg) {
		Population ans;

		std::sort(arg.begin(), arg.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness > rhs.fitness;
		});
		ans = arg;

		ans.resize(2);
		return ans;
	}

	Population ranking(Population arg) {
		Population ans;
		int sum, num;

		std::sort(arg.begin(), arg.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness < rhs.fitness;
		});

		sum = 0;
		for (int i = 0; i < (int)arg.size(); i++) {
			sum = sum + i + 1;
		}
		while (ans.size() < 2) {
			num = (rand() % sum) + 1;
			for (int i = 0; i < (int)arg.size(); i++) {
				num = num - i - 1;
				if (num <= 0) {
					ans.push_back(arg[i]);
					break;
				}
			}
		}

		return ans;
	}

	const funcPtrVec_t call = { temp, ranking };
}

const funcPtr_t selection = Selection::call[SELECTION];

namespace Crossover {
	enum FuncID_e { TEMP };

	Population temp(Population arg) {

		for (unsigned int i = 0; i < arg[i].genome.size(); i++) {
			if (rand() % 2 == 1) {

			}
			else {

				std::swap(arg[0].genome[i], arg[1].genome[i]);
			}
		}

		return arg;
	}

	Population useTable(Population arg) {
		return arg;
	}

	const funcPtrVec_t call = { temp, useTable };
}

const funcPtr_t crossover = Crossover::call[CROSSOVER];

namespace Mutation {
	Population temp(Population arg) {
		for (unsigned int i = 0; i < arg.size(); i++) {
			if (rand() % 10000 < MUTATION_RATE * 10000) {
				for (unsigned int j = 0; j < arg[i].genome.size(); j++) {
					if (rand() % 2 == 1) {
						arg[i].genome[j] = !arg[i].genome[j];
					}
				}
			}
		}
		return arg;
	}

	const funcPtrVec_t call = { temp };
}

const funcPtr_t mutation = Mutation::call[MUTATION];

int main() {
	srand((unsigned)time(NULL));
	Population population(POPULATION_SIZE, GENOME_SIZE);

	population = evaluation(population);
	std::cout << std::endl;
	std::cout << "generation " << 0 << " :" << std::endl;
	population.show();

	for (unsigned int k = 0; k < GENERATION_NUM; k++) {
		Population nextPopulation;

		nextPopulation = population;
		std::sort(nextPopulation.begin(), nextPopulation.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness > rhs.fitness;
		});
		nextPopulation.resize(ELITE_POPULATION_SIZE);
		while (nextPopulation.size() < POPULATION_SIZE) {
			Population temp;

			temp = selection(population);
			temp = crossover(temp);
			temp = mutation(temp);

			std::copy(temp.begin(), temp.end(), back_inserter(nextPopulation));
		}
		nextPopulation.resize(POPULATION_SIZE);

		population = evaluation(nextPopulation);
		std::cout << std::endl;
		std::cout << "generation " << k + 1 << " :" << std::endl;
		population.show();
	}
}