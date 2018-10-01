#include <vector>
#include <iostream>
#include <algorithm>

const int GENOME_SIZE = 32;
const int POPULATION_SIZE = 20;
const int SELECTED_POPULATION_SIZE = 2;
const int GENERATION_NUM = 5;
const double MUTATION_RATE = 0.05;

typedef std::vector<bool> GENOME_t;
typedef int FITNESS_t;

class Individual {
public:
	GENOME_t genome;
	FITNESS_t fitness;
	Individual() {
		for (int i = 0; i < GENOME_SIZE; i++) {
			genome.insert(genome.begin(), rand() % 2);
		}
	}
};

typedef std::vector<Individual> POPULATION_t;

class Population : public POPULATION_t {
public:
	//POPULATION_t population;
	Population() {
		for (int i = 0; i < POPULATION_SIZE; i++) {
			this->push_back(Individual());
		}
	}
	void show() {
		for (POPULATION_t::iterator it = POPULATION_t::begin(); it != POPULATION_t::end(); it++) {
			for (int j = 0; j < GENOME_SIZE; j++) {
				std::cout << it->genome[j];
			}
			std::cout << " " << it->fitness << std::endl;
		}
	}
};

Population evaluation(Population arg) {
	FITNESS_t count;
	for (long unsigned int i = 0; i < arg.size(); i++) {
		count = 0;
		for (int j = 0; j < GENOME_SIZE; j++) {
			if (arg[i].genome[j] == true) count++;
		}
		//std::cout << count << std::endl;
		arg[i].fitness = count;
	}
	return arg;
}

Population selection(Population arg) {
	Population ans;
	/*
	FITNESS_t fitnessSum;
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

Population crossover(Population arg) {
	//Population ans;
	//Individual a, b;

	//ans.resize(0);
	//std::cout << ans.population.size() << std::endl;
	//while (ans.size() < POPULATION_SIZE) {
	for (int i = 0; i < GENOME_SIZE; i++) {
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

Population mutation(Population arg) {
	for (long unsigned int i = 0; i < arg.size(); i++) {
		if (rand() % 10000 < MUTATION_RATE * 10000) {
			std::cout << "mutation!!" << std::endl;
			for (int j = 0; j < GENOME_SIZE; j++) {
				if (rand() % 2 == 1) {
					arg[i].genome[j] = !arg[i].genome[j];
				}
			}
		}
	}
	return arg;
}

Population loop(Population arg) {
	Population b, c, d, e, f;
	/*
	std::cout << std::endl;
	std::cout << "initialized:" << std::endl;
	for(int i = 0; i < (int)a.population.size(); i++){
	for(int j = 0; j < GENOME_SIZE; j++){
	std::cout << a.population[i].genome[j];
	}
	std::cout << " " << a.population[i].fitness << std::endl;
	}
	*/
	b = evaluation(arg);
	std::cout << std::endl;
	std::cout << "evaluated:" << std::endl;
	b.show();
	/*
	for(int i = 0; i < (int)b.population.size(); i++){
	for(int j = 0; j < GENOME_SIZE; j++){
	std::cout << b.population[i].genome[j];
	}
	std::cout << " " << b.population[i].fitness << std::endl;
	}
	*/
	f.resize(0);
	std::cout << std::endl;
	std::cout << "populationSize = " << f.size() << std::endl;
	while (f.size() < POPULATION_SIZE) {
		c = selection(b);
		std::cout << std::endl;
		std::cout << "selecetd:" << std::endl;
		c.show();
		/*
		for(int i = 0; i < (int)c.population.size(); i++){
		for(int j = 0; j < GENOME_SIZE; j++){
		std::cout << c.population[i].genome[j];
		}
		std::cout << " " << c.population[i].fitness << std::endl;
		}
		*/
		d = crossover(c);
		std::cout << std::endl;
		std::cout << "crossovered:" << std::endl;
		d.show();
		/*
		for(int i = 0; i < (int)d.population.size(); i++){
		for(int j = 0; j < GENOME_SIZE; j++){
		std::cout << d.population[i].genome[j];
		}
		std::cout << " " << d.population[i].fitness << std::endl;
		}
		*/
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

int main() {
	Population a;

	for (int k = 0; k < GENERATION_NUM; k++) {
		std::cout << std::endl;
		std::cout << "generation " << k << " :" << std::endl;
		a = loop(a);
	}
}