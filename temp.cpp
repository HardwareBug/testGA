#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>
#include <optional>

std::random_device RANDOM_DEVICE;
std::mt19937 MT(RANDOM_DEVICE());
std::uniform_real_distribution<double> RANDOM_RATE(0.0000, 1.0000);
std::uniform_int_distribution<int> RANDOM_BOOL(0, 1);

const unsigned int GENOME_SIZE = 100;
const unsigned int POPULATION_SIZE = 20;
const unsigned int GENERATION_NUM = 5;
const unsigned int ELITE_POPULATION_SIZE = 2;
const double MUTATION_RATE = 0.01;

// 0 : Two-point Crossover, 1 : Uniform Crossover
const unsigned int TABLE = 1;

const double TABLE_INIT_COST = 1;

// 0 : Temporary, 1 : LifeGame (my function)
const unsigned int EVALUATION = 1;
// 0 : Temporary, 1 : ranking
const unsigned int SELECTION = 1;
// 0 : Temporary, 1 : Use crossover table
const unsigned int CROSSOVER = 1;
// 0 : Temporary
const unsigned int MUTATION = 0;

class Connection {
public:
	unsigned int node1, node2;
	double cost;
	Connection() {};
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
	std::optional<fitness_t> fitness;
	std::optional<fitness_t> pre_fitness;
	connections_t table, choice, cut;
	Individual() {}
	Individual(unsigned int genomeSize) {
		for (unsigned int i = 0; i < genomeSize; i++) {
			genome.insert(genome.begin(), RANDOM_BOOL(MT));

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
			std::cout << " " << it->fitness.value() << std::endl;
		}
	}
};

typedef Population(*funcPtr_t)(Population);
typedef std::vector<funcPtr_t> funcPtrVec_t;

namespace Evaluation {
	Population temp(Population arg) {
		fitness_t count;
		for (unsigned int i = 0; i < arg.size(); i++) {
			arg[i].pre_fitness = arg[i].fitness;
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
			arg[i].pre_fitness = arg[i].fitness;
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
						arg[i].fitness.value() += 1;
					}
					else {
						arg[i].fitness.value() -= 1;
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

		while (ans.size() < 2) {
			sum = 0;
			for (int i = 0; i < (int)arg.size(); i++) {
				sum = sum + i + 1;
			}
			std::uniform_int_distribution<int> random_temp(1, sum);
			num = random_temp(MT);
			for (int i = 0; i < (int)arg.size(); i++) {
				num = num - i - 1;
				if (num <= 0) {
					ans.push_back(arg[i]);
					arg.erase(arg.begin() + i);
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

		for (unsigned int i = 0; i < arg[0].genome.size(); i++) {
			if (RANDOM_BOOL(MT) == 1) {

			}
			else {

				std::swap(arg[0].genome[i], arg[1].genome[i]);
			}
		}

		return arg;
	}

	Population useTable(Population arg) {
		connections_t table, choice, cut;
		std::copy(arg[0].table.begin(), arg[0].table.end(), back_inserter(table));
		
		std::vector<std::vector<unsigned int> > nodes;
		for(unsigned int i = 0; i < arg[0].genome.size(); i++){
			std::vector<unsigned int> temp;
			temp.push_back(i);
			nodes.push_back(temp);
		}

		//
		//todo (change graph partitioning algorithm)
		//

		while(nodes.size() > 2){
			double sum, randNum;
			sum = 0;
			for (int i = 0; i < (int)table.size(); i++) {
				sum += table[i].cost;
			}
			if(sum != 0){
				std::uniform_real_distribution<double> random_real_temp(0, sum);
				randNum = random_real_temp(MT);
				//std::cout << randNum << std::endl;
				for (int i = 0; i < (int)table.size(); i++) {
					randNum -= table[i].cost;
					if(table[i].cost != 0){
						if(randNum <= 0){

							unsigned int node1_p, node2_p;
							for(unsigned int j = 0; j < nodes.size(); j++){
								for(unsigned int k = 0; k < nodes[j].size(); k++){
									if( nodes[j][k] == table[i].node1){
										node1_p = j;
									}
									if( nodes[j][k] == table[i].node2){
										node2_p = j;
									}
								}
							}
							if(node1_p != node2_p){
								if(node1_p > node2_p){
									std::swap(node1_p, node2_p);
								}
								std::copy(nodes[node2_p].begin(), nodes[node2_p].end(), back_inserter(nodes[node1_p]));
								nodes[node2_p].clear();
								nodes.erase(nodes.begin() + node2_p);
							}

							choice.push_back(table[i]);
							table.erase(table.begin() + i);
							break;
						}
					}
				}
			}else{
				std::uniform_int_distribution<int> random_int_temp(0, table.size()-1);
				randNum = random_int_temp(MT);

				unsigned int node1_p, node2_p;
				for(unsigned int j = 0; j < nodes.size(); j++){
					for(unsigned int k = 0; k < nodes[j].size(); k++){
						if( nodes[j][k] == table[randNum].node1){
							node1_p = j;
						}
						if( nodes[j][k] == table[randNum].node2){
							node2_p = j;
						}
					}
				}
				if(node1_p != node2_p){
					if(node1_p > node2_p){
						std::swap(node1_p, node2_p);
					}
					std::copy(nodes[node2_p].begin(), nodes[node2_p].end(), back_inserter(nodes[node1_p]));
					nodes[node2_p].clear();
					nodes.erase(nodes.begin() + node2_p);
				}

				choice.push_back(table[randNum]);
				table.erase(table.begin() + randNum);
			}
			/*
			for (unsigned int i = 0; i < nodes.size(); i++) {
				for (unsigned int j = 0; j < nodes[i].size(); j++) {
					std::cout << nodes[i][j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
			*/
		}
		
		for (unsigned int i = 0; i < nodes[0].size(); i++) {
			std::swap(arg[0].genome[nodes[0][i]], arg[1].genome[nodes[0][i]]);
		}
		
		for (unsigned int i = 0; i < nodes.size(); i++) {
			for (unsigned int j = 0; j < nodes[i].size(); j++) {
				std::cout << nodes[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		return arg;
	}

	const funcPtrVec_t call = { temp, useTable };
}

const funcPtr_t crossover = Crossover::call[CROSSOVER];

namespace Mutation {
	Population temp(Population arg) {
		for (unsigned int i = 0; i < arg.size(); i++) {
			if (RANDOM_RATE(MT) < MUTATION_RATE) {
				for (unsigned int j = 0; j < arg[i].genome.size(); j++) {
					if (RANDOM_BOOL(MT) == 1) {
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
	//srand((unsigned)time(NULL));
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
		
		//
		//todo (add UpdateCrossoverTable)
		//
		
		std::cout << std::endl;
		std::cout << "generation " << k + 1 << " :" << std::endl;
		population.show();
	}
}