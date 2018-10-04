#include <vector>
#include <iostream>
#include <algorithm>

const unsigned int GENOME_SIZE = 32;
const unsigned int POPULATION_SIZE = 20;
const unsigned int GENERATION_NUM = 5;
const double MUTATION_RATE = 0.05;

const unsigned int SELECTED_POPULATION_SIZE = 2;

typedef std::vector<bool> genome_t;
typedef int fitness_t;

class Individual {
public:
	genome_t genome;
	fitness_t fitness;
	Individual() {}
	Individual(unsigned int genomeSize) {
		for (unsigned int i = 0; i < genomeSize; i++) {
			genome.insert(genome.begin(), rand() % 2);
		}
	}
};

typedef std::vector<Individual> population_t;

class Population : public population_t {
public:
	//population_t population;
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
}

namespace Selection {
	Population temp(Population arg) {
		Population ans;
		/*
		fitness_t fitnessSum;
		std::vector<double> probabilityAccumulation(POPULATION_SIZE)
		double randNum;
		for(int i = 0; i < POPULATION_SIZE; i++){
		fitnessSum += arg.population[i].fitness;
		probabilityAccumulation[i] = arg.population[i].fitness;
		}
		for(int i = 0; i < POPULATION_SIZE; i++){
		probabilityAccumulation[i] /= fitnessSum;
		}
		while(ans.size < SELECTED_POPULATION){
		randNum = rand() % 1;
		for(int i = 0; i < arg.size(); i++){
		if(randNum < probabilityAccumulation){
		ans.push_back(arg.population[i]);
		}
		}
		}
		*/
		std::sort(arg.begin(), arg.end(), [](const Individual &lhs, const Individual &rhs) {
			return lhs.fitness > rhs.fitness;
		});
		ans = arg;
		ans.resize(SELECTED_POPULATION_SIZE);
		return ans;
	}
}

namespace Crossover {
	Population temp(Population arg) {
		//Population ans;
		//Individual a, b;

		//ans.resize(0);
		//std::cout << ans.population.size() << std::endl;
		//while (ans.size() < POPULATION_SIZE) {
		for (unsigned int i = 0; i < arg[i].genome.size(); i++) {
			if (rand() % 2 == 1) {
				//std::cout << "1" << std::endl;
				//a.genome[i] = arg[0].genome[i];
				//b.genome[i] = arg[1].genome[i];
			}
			else {
				//std::cout << "0" << std::endl;
				//a.genome[i] = arg[1].genome[i];
				//b.genome[i] = arg[0].genome[i];
				std::swap(arg[0].genome[i], arg[1].genome[i]);
			}
		}
		//ans.push_back(a);
		//ans.push_back(b);
		//}
		return arg;
	}
}

namespace Mutation {
	Population temp(Population arg) {
		for (unsigned int i = 0; i < arg.size(); i++) {
			if (rand() % 10000 < MUTATION_RATE * 10000) {
				std::cout << "mutation!!" << std::endl;
				for (unsigned int j = 0; j < arg[i].genome.size(); j++) {
					if (rand() % 2 == 1) {
						arg[i].genome[j] = !arg[i].genome[j];
					}
				}
			}
		}
		return arg;
	}
}

/*
Population loop(Population arg) {
Population b, c, d, e, f;
std::cout << std::endl;
std::cout << "initialized:" << std::endl;
for(int i = 0; i < (int)a.population.size(); i++){
for(int j = 0; j < GENOME_SIZE; j++){
std::cout << a.population[i].genome[j];
}
std::cout << " " << a.population[i].fitness << std::endl;
}
b = evaluation(arg);
std::cout << std::endl;
std::cout << "evaluated:" << std::endl;
b.show();
for(int i = 0; i < (int)b.population.size(); i++){
for(int j = 0; j < GENOME_SIZE; j++){
std::cout << b.population[i].genome[j];
}
std::cout << " " << b.population[i].fitness << std::endl;
}
//f.resize(0);
std::cout << std::endl;
std::cout << "populationSize = " << f.size() << std::endl;
while (f.size() < b.size()) {
c = selection(b);
std::cout << std::endl;
std::cout << "selecetd:" << std::endl;
c.show();
for(int i = 0; i < (int)c.population.size(); i++){
for(int j = 0; j < GENOME_SIZE; j++){
std::cout << c.population[i].genome[j];
}
std::cout << " " << c.population[i].fitness << std::endl;
}
d = crossover(c);
std::cout << std::endl;
std::cout << "crossovered:" << std::endl;
d.show();
for(int i = 0; i < (int)d.population.size(); i++){
for(int j = 0; j < GENOME_SIZE; j++){
std::cout << d.population[i].genome[j];
}
std::cout << " " << d.population[i].fitness << std::endl;
}
e = mutation(d);
std::cout << std::endl;
std::cout << "mutated:" << std::endl;
e.show();

copy(e.begin(), e.end(), back_inserter(f));
std::cout << std::endl;
std::cout << "populationSize = " << f.size() << std::endl;
}
return f;
}
*/

typedef Population(*gaFuncPtr_t)(Population);
typedef std::vector<gaFuncPtr_t> gaFuncPtrVec_t;

class MyGA {
private:
	unsigned int genomeSize = GENOME_SIZE;
	unsigned int populationSize = POPULATION_SIZE;
	unsigned int generationNum = GENERATION_NUM;
	double mutationRate = MUTATION_RATE;

	int evalutaionFuncID = 0;
	int selectionFuncID = 0;
	int crossoverFuncID = 0;
	int mutationFuncID = 0;

	gaFuncPtrVec_t evaluationFuncPtrVec = { Evaluation::temp };
	gaFuncPtrVec_t selectionFuncPtrVec = { Selection::temp };
	gaFuncPtrVec_t crossoverFuncPtrVec = { Crossover::temp };
	gaFuncPtrVec_t mutationFuncPtrVec = { Mutation::temp };

	Population evaluation(Population arg) {
		return evaluationFuncPtrVec[evalutaionFuncID](arg);
	}

	Population selection(Population arg) {
		return selectionFuncPtrVec[selectionFuncID](arg);
	}

	Population crossover(Population arg) {
		return crossoverFuncPtrVec[crossoverFuncID](arg);
	}

	Population mutation(Population arg) {
		return mutationFuncPtrVec[mutationFuncID](arg);
	}
public:
	void done() {
		Population population(POPULATION_SIZE, GENOME_SIZE), populationNext, temp;

		population = this->evaluation(population);
		std::cout << std::endl;
		std::cout << "generation " << 0 << " :" << std::endl;
		population.show();

		for (unsigned int k = 0; k < GENERATION_NUM; k++) {
			populationNext.resize(0);
			while (populationNext.size() < populationSize) {
				temp = this->selection(population);
				temp = this->crossover(temp);
				temp = this->mutation(temp);

				copy(temp.begin(), temp.end(), back_inserter(populationNext));
			}
			population = this->evaluation(population);
			std::cout << std::endl;
			std::cout << "generation " << k + 1 << " :" << std::endl;
			population.show();
		}
	}

	void setGenomeSize(unsigned int size) {
		genomeSize = size;
	}

	void setPopulationSize(unsigned int size) {
		populationSize = size;
	}

	void setGenerationNum(unsigned int num) {
		generationNum = num;
	}

	void setMutationRate(double rate) {
		mutationRate = rate;
	}
};

int main() {
	//Population a(POPULATION_SIZE, GENOME_SIZE);
	MyGA ga;
	ga.done();
	/*
	for (unsigned int k = 0; k < GENERATION_NUM; k++) {
	std::cout << std::endl;
	std::cout << "generation " << k << " :" << std::endl;
	a = loop(a);
	}
	*/
}