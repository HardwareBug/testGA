#define EIGEN_NO_DEBUG
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_MPL2_ONLY
#include "Eigen/Dense"

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

typedef std::vector<bool> genome_t;
typedef int fitness_t;

class Individual {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	genome_t genome;
	std::optional<fitness_t> fitness;
	std::optional<fitness_t> pre_fitness;
	Eigen::MatrixXd table;
	
	Individual() {}
	Individual(unsigned int genomeSize) {
		table.resize(GENOME_SIZE, GENOME_SIZE);
		for (int i = 0; i < genomeSize; i++) {
			genome.insert(genome.begin(), RANDOM_BOOL(MT));

			for (int j = 0; j < genomeSize; j++) {
				if (TABLE == 0) {
					if (std::abs(i - j) == 1 || std::abs(i - j) == GENOME_SIZE - 1) {
						table(i, j) = TABLE_INIT_COST;
					}else{
						table(i, j) = 0;
					}
				}
				else if (TABLE == 1) {
					if (std::abs(i - j) != 0) {
						table(i, j) = TABLE_INIT_COST;
					}else{
						table(i, j) = 0;
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
		
		Eigen::MatrixXd degree(GENOME_SIZE, GENOME_SIZE);
		degree = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		for (int i = 0; i < GENOME_SIZE; i++) {
			double temp = 0;
			
			for (int j = 0; j < GENOME_SIZE; j++) {
				temp += arg[0].table(i, j);
			}
			
			degree(i, i) = temp;
		}
		//std::cout << degree << std::endl << std::endl;
		
		//
		//
		//
		
		Eigen::MatrixXd laplacian(GENOME_SIZE, GENOME_SIZE);
		laplacian = degree - arg[0].table;
		//std::cout << laplacian << std::endl << std::endl;
		
		//
		//
		//
		
		Eigen::MatrixXd degreeLeft(GENOME_SIZE, GENOME_SIZE);
		degreeLeft = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		for (int i = 0; i < GENOME_SIZE; i++) {
			degreeLeft(i, i) = 1 / std::sqrt(degree(i, i));
		}
		//std::cout << degreeLeft << std::endl << std::endl;
		
		Eigen::MatrixXd degreeRight(GENOME_SIZE, GENOME_SIZE);
		degreeRight = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		for (int i = 0; i < GENOME_SIZE; i++) {
			degreeRight(i, i) = std::sqrt(degree(i, i));
		}
		//std::cout << degreeRight << std::endl << std::endl;
		
		Eigen::MatrixXd target(GENOME_SIZE, GENOME_SIZE);
		target = degreeLeft * laplacian * degreeRight;
		//std::cout << target << std::endl << std::endl;
		
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(target);
		Eigen::MatrixXd value, vector;
		value = eig.eigenvalues();
		vector = eig.eigenvectors();
		//std::cout << value << std::endl << std::endl;
		//std::cout << vector << std::endl << std::endl;
		
		for (unsigned int i = 0; i < GENOME_SIZE; i++) {
			if(vector(1, i) < 0){
				std::swap(arg[0].genome[i], arg[1].genome[i]);
				std::cout << i << " ";
			}
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