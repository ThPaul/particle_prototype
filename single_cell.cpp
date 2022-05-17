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

std::vector<Utils::Vector3d> one_with_all_after(const ParticleVector& ParVec ,const int& Startindex){
	std::vector<Utils::Vector3d> forces(ParVec.size());
	auto& FirstElem = ParVec[Startindex];
	for(int i=Startindex+1;i<ParVec.size();++i){
		auto force = FirstElem.pos()-ParVec[i].pos();
		forces[i] += force;
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
	std::vector<hpx::shared_future<std::vector<Utils::Vector3d>>> futures;
	futures.reserve(n);
	for (int i=0; i<n ; ++i){
		particleList[i].pos()={i,0,0};	}	
	for (int i=0; i<n ;++i){

		futures.push_back(hpx::async(one_with_all_after<typeof(particleList)>,particleList,i));	

	}
	for (int i=0; i<n; ++i){
		Utils::Vector3d sum;
		sum.fill(0);
		for (int j=0; j<i; ++j){
			sum-=futures[j].get()[i];
		}
		for (int j=i+1;j<n;++j){
			sum+=futures[i].get()[j];
		}
		particleList[i].force()=sum;
		if(n==50){
			std::cout<<i<<" "<<sum[0]<<std::endl; 
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
