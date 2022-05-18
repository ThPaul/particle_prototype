#include <stdlib.h>

#include "config.hpp"
#include "all_particles.hpp"
#include "Cell.hpp"
#include "utils/for_each_pair.hpp"
#include "utils/Vector.hpp"
#include "velocity_verlet_inline.hpp"
#include "timer.hpp"

#include <iostream>
#include <string>
#include <hpx/hpx.hpp>
#include <hpx/include/future.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/actions.hpp>
#include <vector>


template <typename Particle>
Cell<Particle> setup_cell(int n) {
	Cell<Particle> cell;
	for (int i=0;i<n;i++) cell.particles().insert(Particle{});
	return cell; 
};


template <typename ParticleVector>

std::vector<Utils::Vector3d> one_with_all_after(const ParticleVector& ParVec ,const int& Startindex, const int& Endindex){
	std::vector<Utils::Vector3d> forces(ParVec.size());

	for(int i=Startindex; i<= Endindex; ++i){
		// possible to optimize by decreasing size of vector and obmmitting Zeros
		auto& FirstElemPos = ParVec[i].pos();
		for(int j=i+1;j<ParVec.size();++j){
			auto force = FirstElemPos-ParVec[j].pos();
			forces[i]+=force;
			forces[j] -= force;
		}

	}
	return forces;

}

template <typename ParticleIterable>

void all_with_all(ParticleIterable& particles) {
	int n=particles.size();
	std::vector<typeof(*particles.begin())> particleList;
	for (auto it = particles.begin(); it!= particles.end(); ++it){
		particleList.push_back(*it);	
	}	
	for (int i=0; i<n ; ++i){
		particleList[i].pos()={double(i),0,0};	}	
	int NrOfPar = 20;
	if(NrOfPar>n) NrOfPar=n;
	int parSize = n / NrOfPar;
	int NrOfLarger = n % NrOfPar ;
	auto Index1d = [parSize, NrOfLarger](int par, int i) { return par*parSize + i + ((par<NrOfLarger) ? par : NrOfLarger );};
/* NOT USED
	auto Index2d = [parSize, NrOfLarger](int i) {
		if (i/(parSize+1) < NrOfLarger){
			std::array result{i/(parSize+1),i%(parSize+1)};
			return result;
		}
		else{
			int r = i - (parSize+1)*NrOfLarger;
			std::array result{r/parSize+NrOfLarger,r%parSize};
			return result;
		}
	}; */

	/*for (int i=0; i<20; ++i){
	  int random=rand()%n;
	  auto b=Index2d(random);
	  auto c=Index1d(b[0],b[1]);

	  std::cout << n << "NUMBER " << random  << " 2d " << b[0] <<"," << b[1] << " back " << c <<  std::endl;

	  }*/	

	std::vector<hpx::shared_future<std::vector<Utils::Vector3d>>> allPartitionsf;
	allPartitionsf.reserve(NrOfPar);
	for (int i=0; i<NrOfPar ; ++i){
		int start = Index1d(i,0);
		int end = Index1d(i+1,0)-1;	
		allPartitionsf.push_back(hpx::async(one_with_all_after<typeof(particleList)>,particleList,start,end));	
	}

	std::vector<Utils::Vector3d> finalForces(n);
	std::vector<Utils::Vector3d> tempForces(n);
	for (int i=0; i<allPartitionsf.size(); ++i){
		tempForces=allPartitionsf[i].get();
		for (int j=0; j<n; ++j){
		finalForces[j]+=tempForces[j];
		}
	}
	  for (int i=0; i<n; ++i){
	   particleList[i].force()=finalForces[i];
	  if(n==50){
	  std::cout<<i<<" "<<finalForces[i][0]<<std::endl; 
	  }
	  }

}


template <typename ParticleIterable> 
void single_pass(ParticleIterable& particles) {
	velocity_verlet_step_1(particles, 0.01);
}


template <typename Particle>
void time(int n) {
	auto cell = setup_cell<Particle>(n);
	auto tick = Timer::now();
	single_pass(cell.particles());
	auto single_pass_time = Duration(Timer::now() -tick);

	tick = Timer::now();
	all_with_all(cell.particles());
	auto all_with_all_time = Duration(Timer::now() -tick);
	//  std::cerr << (cell.particles().begin())->pos() <<std::endl;
	std::cout <<n <<" " <<to_ms(single_pass_time) << " "<<to_ms(all_with_all_time)<<std::endl;
}


template <typename Particle>
void run_timing(std::string particle_name) {
	std::cout << particle_name <<" Size: "<<sizeof(Particle)<< std::endl;
	for (int n: {50, 5000, 10000, 15000, 20000,25000,30000,35000,40000, 45000, 50000, 55000, 60000, 65000, 70000, 75000, 80000}) {
		time<Particle>(n);
	}
}
int main(int argc, char** argv) {
	//  run_timing<CurrentEsParticle>("Current Espresso Particle");
	//  run_timing<MinimalFlatParticle<456>>("Re-ordered properties");
	run_timing<MinimalFlatParticle<0>>("Minimal flat particle");
	//  run_timing<SoABackedParticle>("Struct-of-arrays backed particle");
	return 0;
}
