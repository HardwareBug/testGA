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

std::uniform_real_distribution<double> RANDOM_TABLE(0.0, 0.1);

const unsigned int GENOME_SIZE = 100;
const unsigned int POPULATION_SIZE = 20;
const unsigned int GENERATION_NUM = 100;
const unsigned int ELITE_POPULATION_SIZE = 2;
const double MUTATION_RATE = 0.01;

const unsigned int UPDATE_TABLE_DATASIZE_LINE = 100;
const unsigned int UPPER_DATASIZE = 20;
const unsigned int DOWNER_DATASIZE = 20;
const unsigned int USE_UPDATE_TABLE = 1;

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
	Eigen::MatrixXd connect;
	Eigen::MatrixXd cut;
	bool has_data;
	
	Individual() {}
	Individual(unsigned int genomeSize) {
		table.resize(GENOME_SIZE, GENOME_SIZE);
		connect.resize(GENOME_SIZE, GENOME_SIZE);
		cut.resize(GENOME_SIZE, GENOME_SIZE);
		has_data = false;
		/*
		table(0, 0) = 0.0;
		table(0, 1) = 0.8;
		table(0, 2) = 0.1;
		table(0, 3) = 0.1;
		table(0, 4) = 0.1;
		table(1, 0) = 0.8;
		table(1, 1) = 0.0;
		table(1, 2) = 0.4;
		table(1, 3) = 0.2;
		table(1, 4) = 0.2;
		table(2, 0) = 0.1;
		table(2, 1) = 0.4;
		table(2, 2) = 0.0;
		table(2, 3) = 0.9;
		table(2, 4) = 0.7;
		table(3, 0) = 0.1;
		table(3, 1) = 0.2;
		table(3, 2) = 0.9;
		table(3, 3) = 0.0;
		table(3, 4) = 0.8;
		table(4, 0) = 0.1;
		table(4, 1) = 0.2;
		table(4, 2) = 0.7;
		table(4, 3) = 0.8;
		table(4, 4) = 0.0;
		for (int i = 0; i < genomeSize; i++) {
			genome.insert(genome.begin(), RANDOM_BOOL(MT));
		}
		*/
		
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
			if(it->pre_fitness.has_value()){
				std::cout << " " << it->fitness.value() << " " << it->pre_fitness.value() << std::endl;
			}else{
				std::cout << " " << it->fitness.value() << std::endl;
			}
		}
	}
};

typedef Population(*funcPtr_t)(Population);
typedef std::vector<funcPtr_t> funcPtrVec_t;

namespace Evaluation {
	Population temp(Population arg) {
		fitness_t count;
		for (unsigned int i = 0; i < arg.size(); i++) {
			arg[i].pre_fitness.swap(arg[i].fitness);
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
			arg[i].pre_fitness.swap(arg[i].fitness);
			
			arg[i].fitness = 0;
			for (int j = 0; j < (int)arg[i].genome.size(); j++) {
				count = 0;
				for (int k = -1; k <= 1; k++) {
					for (int l = -1; l <= 1; l++) {
						if (0 <= j + k + l*size && j + k + l*size < (int)arg[i].genome.size() && 0 <= j%size + k && j%size + k < size && 0 <= j/size + l && j/size + l < size) {
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
		Eigen::MatrixXd adjacency(GENOME_SIZE, GENOME_SIZE);
		adjacency = arg[0].table;
		for (int i = 0; i < GENOME_SIZE; i++) {
			for (int j = 0; j < GENOME_SIZE; j++) {
				if(i != j){
					adjacency(i, j) += RANDOM_TABLE(MT);
				}
			}
		}
		//std::cout << adjacency << std::endl << std::endl;
		
		Eigen::MatrixXd degree(GENOME_SIZE, GENOME_SIZE);
		degree = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		for (int i = 0; i < GENOME_SIZE; i++) {
			double temp = 0;
			
			for (int j = 0; j < GENOME_SIZE; j++) {
				temp += adjacency(i, j);
			}
			
			degree(i, i) = temp;
		}
		//std::cout << degree << std::endl << std::endl;
		
		//
		//
		//
		
		Eigen::MatrixXd laplacian(GENOME_SIZE, GENOME_SIZE);
		laplacian = degree - adjacency;
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
			degreeRight(i, i) = 1 / std::sqrt(degree(i, i));
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
			if(vector(i, 1) < 0){
				std::swap(arg[0].genome[i], arg[1].genome[i]);
				//std::cout << i << " ";
			}
		}
		//std::cout << std::endl;
		
		for (unsigned int i = 0; i < GENOME_SIZE; i++) {
			for (unsigned int j = 0; j < GENOME_SIZE; j++) {
				if(i != j && (vector(i, 1) < 0) && (vector(j, 1) < 0)){
					arg[0].connect(i, j) = 1;
					arg[1].connect(i, j) = 1;
				}else{
					arg[0].connect(i, j) = 0;
					arg[1].connect(i, j) = 0;
				}
			}
		}
		arg[0].cut = Eigen::MatrixXd::Ones(GENOME_SIZE, GENOME_SIZE) - Eigen::MatrixXd::Identity(GENOME_SIZE, GENOME_SIZE) - arg[0].connect ;
		arg[1].cut = Eigen::MatrixXd::Ones(GENOME_SIZE, GENOME_SIZE) - Eigen::MatrixXd::Identity(GENOME_SIZE, GENOME_SIZE) - arg[1].connect ;
		
		//std::cout << arg[0].connect << std::endl << std::endl;
		//td::cout << arg[0].cut << std::endl << std::endl;
		
		arg[0].has_data = true;
		arg[1].has_data = true;
		
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

class CrossoverTable {
private:
	Population data;
	
	Population update_table(Population arg){
		std::sort(data.begin(), data.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness > rhs.fitness;
		});
		
		Eigen::MatrixXd upper_connect = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd upper_cut = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd downer_connect = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd downer_cut = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd connect = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd cut = Eigen::MatrixXd::Zero(GENOME_SIZE, GENOME_SIZE);
		Eigen::MatrixXd teble = arg[0].table;
		
		for(unsigned int i = 0; i < UPPER_DATASIZE; i++){
			upper_connect = upper_connect + data[i].connect;
			upper_cut = upper_cut + data[i].cut;
		}
		upper_connect = upper_connect.array() / UPPER_DATASIZE;
		upper_cut = upper_cut.array() / UPPER_DATASIZE;
		//std::cout << upper_connect << std::endl << std::endl;
		//std::cout << upper_cut << std::endl << std::endl;
		
		for(unsigned int i = data.size() - DOWNER_DATASIZE; i < data.size(); i++){
			downer_connect = downer_connect + data[i].connect;
			downer_cut = downer_cut + data[i].cut;
		}
		downer_connect = downer_connect.array() / DOWNER_DATASIZE;
		downer_cut = downer_cut.array() / DOWNER_DATASIZE;
		//std::cout << downer_connect << std::endl << std::endl;
		//std::cout << downer_cut << std::endl << std::endl;
		
		connect = upper_connect - downer_connect;
		cut = upper_cut - downer_cut;
		//std::cout << connect << std::endl << std::endl;
		//std::cout << cut << std::endl << std::endl;
		
		connect = connect.array() + 1;
		cut = cut.array() * (-1) + 1;
		
		teble = arg[0].table.array() * connect.array() * cut.array();
		if(USE_UPDATE_TABLE){
			for(unsigned int i = 0; i < arg.size(); i++){
				arg[i].table = teble;
			}
		}
		
		data.erase(data.begin(), data.end());
		
		return arg;
	}
public:
	Population update(Population arg){
		Population temp = arg;
		
		for(unsigned int i = 0; i < arg.size(); i++){
			if(temp[i].has_data){
				temp[i].has_data = false;
				temp[i].fitness = temp[i].fitness.value() - temp[i].pre_fitness.value();
				temp[i].pre_fitness.reset();
				data.push_back(temp[i]);
			}
		}
		if(data.size() >= UPDATE_TABLE_DATASIZE_LINE){
			arg = update_table(arg);
		}
		
		return arg;
	}
	
	void show(){
		data.show();
	}
};

const funcPtr_t mutation = Mutation::call[MUTATION];

int main() {
	//srand((unsigned)time(NULL));
	Population population(POPULATION_SIZE, GENOME_SIZE);
	CrossoverTable crossoverTable;

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
		
		population = crossoverTable.update(population);
		
		std::cout << std::endl;
		std::cout << "generation " << k + 1 << " :" << std::endl;
		population.show();
	}
	
	std::cout << std::endl;
	std::cout << "crossoverTable :" << std::endl;
	crossoverTable.show();
}