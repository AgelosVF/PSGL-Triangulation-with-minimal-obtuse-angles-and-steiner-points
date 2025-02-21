Created by: ΑΓΓΕΛΟΣ ΒΕΝΕΤΗΣ ΦΑΝΟΥΡΑΚΗΣ sdi1300022
My solutions for the 7th Computational Geometry Challenge (2025)

Files and Directories:

├── build 	<-- Inside will be built the executable opt_triangulation
├── Header Files 	<-- Contrains the header files
│   ├── boost_utils.hpp		<-- Functions used to read from json file	
│   ├── CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h <--Custom cdt class used.Added more methods to avoid
│   │								flips.
│   ├── BoundaryType.hpp	<-- contains function that are used for recognizing the type of triangulation
│   ├── SteinerPoints.hpp	<-- contains the functions that are used across methods for inserting steiner
│   │					points centroid/circumcenter/projection/merge/longest_edge
│   ├── random.hpp	<-- contains functions that use randomness to insert steiner points to the triangulation
│   └── triangulation_utils.hpp		<--Functions rest of functions used between the files
├── make		<-- Inside there is the CMakeLists.txt
│   ├── CMakeLists.txt
│   └── Lib   <-- Empty directory where you can run cmake .. to compile
└── Sources		<-- All the source code for the functions used.
    ├── boost_utils.cpp		<--Source code of functions used to read the json file.
    ├── main.cpp		<--Source code of the main function.
    ├── outputfile.cpp		<--Source code of the function used to create the output.
    ├── SteinerPoints.cpp	<--Source code of the functions used to find and add steiner points.
    ├── triangulation_utils.cpp <--Source code of other helping functions.
    ├── random.cpp		<--Source code of the functions that use randomness to insert steiner.
    ├── BoundaryType.cpp	<--Source code of functions that recognize the category of problem.
    ├── triangulation_utils.cpp <--Source code of other helping functions.
    │			   	   to recreate the triangulation produced by the previous project.
    ├── LocalSearch.cpp		<--Source code of the local search.
    ├── SimulatedAnnealing.cpp	<--Source code of the sumulated annealing.
    └── AntCollony		<--Source code of the ant collony method.
How to compile and run:
	Inside the make directory there is the CMakeLists.txt and another directory called Lib. If you go in the
	lib directory and call cmake .. the makefile will be created in the same directory along with some other
	files. After that you need to call make and the executable opt_triangulation will be created in the 
	build directory.
	To run the program you call "./opt_triangulation -i /path/to/input.json -o /path/to/output.json" and you
	can use the flag -preselected_params in order for it to automaticaly choose the method and parameters.

How boundary type works:
	For convex no constrains we use the cgal is_convex()
	For convex open/closed constrains we put the boundary in a graph. Then we put each additional edge
	and check if starting from it there is a cycle containing it. If there is we got a cycle that may involve
	the boundary but for sure involves an additional constraint edge.
	For non convex with parallel x,y to xx' and yy' we check if all the connected points have either the same
	x or y.
	
How LocalSearch works:
	Local search goes through all the obtuse faces of triangulation and for each test each of the 5
	methods to see how many obtuse faces would be after adding a steiner. All methods use methods of
	the parent class constrained_triangulation in order not to trigger any unwanted flips durring the
	editing of the cdt. Each method has the option to be used to a copy of the cdt so we can test it
	and the actual cdt. Also they return the number of obtuse triangles after the steiner is applied.
	The 5 methods are:
		-steiner on the midpoint of the longest edge
		-steiner on the projection of the obtuse angle
		-steiner on the centroid
		-steiner on the circumcenter:
			First we check if the circumcenter of the face lands inside a neighbor. If it does then
			we check if the common edge is constrained. If it isn't we add the steiner inside the 
			neighbor and then we make the steiner-obtuse_vertex edge a constrain. This way the old
			edge in between is removed and the edges we wanted are formed. After that we remove the
			constrain. For all the adding and removing constrains we use the methods of the parent
			class constrained_triangulation in order no to triger any unwanted delaunay flips.
		-steiner on the centroid of merged obtuse neighbors:
			For this steiner we got 2 function one called simulate_merge and another that adds the
			steiner. 
			The simulate checks with each neighbor and if we can merge them in a copy of the cdt.
			If the merge is allowed it: 
				1. removes the all the edges of the faces we will merge by removing their points.
				2. adds the steiner on the centroid of the polygon the merged faces would make.
				3. adds back all the removed edges and makes the polygon a constrain.
				4. removes the constrains.
				5. counts the obtuse faces after the merge and if that number is the best up
				to here it saves a point to indicate the neighbor that was used for the merge.
				6. resets the triangulation back to its previous state and check the other 
				neighbors.
			That function is called by the actual merge function which if a viable neighbor is returned
			to it, it performs the merge with it.
			

	After checking all the methods we choose the one which had the best result. In case of ties the order
	goes centroid,projection longest edge, circumcenter merge. Then we check if the obtuse count we got
	is less than the one we had before. If it is we apply the method.
	This loops for L times(from input) or until we checked each face and we didnt find any method that would
	reduce the count.

How Simulated Annealing works:
	The methods for steiner points used are almost identical with some technical changes on the returns of 
	the functions and also the circumcenter and merge  methods in case they dont work they do a steiner on
	the centroid.
	The function loops until the Temperature reachs 0 or bellow or L times. In each loop it goes through
	each face and chooses 1 of the methods at random. Then it counts how many obtuse faces we would have
	if we applied that method. We use that number to calculate the energy (test_E=test_obtuse*a +steiner*b 
	and that to calculate the delta energy. If delta energy<0 or a with a chance of e^{-delta_energy/Temp}
	we apply the steiner point. Then if the current enery is less than the best enery we save the current
	cdt as the best.

How Ant colony works:
	 The main function of the ant colony is the ant_colony(). It starts bt creating vectors for the: 
	pheromones,ant_energy and ant_method. After it starts the loop. In the loop for each ant we create
	a a thread which calls ImproveTriangulation(). 
	 In ImproveTriangulation() a random face among the obtuse is selected and a method. For the method we
	use the heuristic functions and the pheromones vector to create a weighted random selection with each
	method having p(k)= (t^xi * h^psi)/ Σ_{i} t^x_{i}h^ps_{i} chance. After the method is selected we 
	create a copy of the cdt and use that method to insert a steiner in the random obtuse face we selected
	before. Each method returns a pointer to the neighboring faces it also changed. Those neighboring faces
	and face that were changed by the steiner insertion are returned as a vector pair {face,neighbor} and
	are saved in a vector called ant_faces.Also the selected_method and energy are saved to the vectors to
	which pointers to are passed as arguments.
	 Back in ant_colony we wait for each thread to finish. Since each thread uses its own copy of the cdt
	and only updates its own vectors we dont need to use a lock since we got no critical section. After 
	each thread (ant) has finished we create a new vector in which we store each energy and its index.
	In the vectors ant_energy, ant_method and ant_faces each ant has its own index for example the second
	ant has its energy stored at ant_energy[1], the method it used at ant_method[1] and the faces it changed
	at ant_faces[1].So after we create the energy-index pair vector we sort it by energy to know which 
	ants to prioritize. After that we call the resolve_conflicts() function to filter the ants.
	 The resolve_conflicts function takes the energy_index_pairs and ant_faces vectors as arguments.
	Then starting from the first ant index in the energy_index_pairs (the one with the best energy) it goes
	throught the other ants and if ant of the faces they changed were the same it changes their index in
	the energy_index_pairs to -1 so they can be skipped. This none of the remaining ants have any common
	changed face and when 2 or more had some common we kept the one with the better energy.
	  After that back at the ant_colony() its time to apply the steiner points of the remaining ants.
	We go throught the energy_index_pairs vector to get the indexes of the ants starting from the one
	with the lower energy and skipping the ones whose index is -1. If the energy is better from the energy
	of the previous cycle we apply the steiner. If the ants energy was 0 it means it reduced the obtuse 
	faces to 0 se we terminate. After that we increase or reduce the methods pheromon if it was applied or
	not. After going throught all the ants we update the current number of steiners and obtuse faces in the
	cdt, update its energy and go back at the start of the loop. The loop goes L times.
	 The functions for the steiner points used are almost identical to the ones used in simulated annealing.
	The other functions are the heuristics for each method and the increase/decrease pheromones which were 
	taken from the class notes.
How reduce_random_local works:
	  The reduce_random_local is a mix of local search with simulated annealing. It starts by creating a
	 copy of the cdt and adds a steiner point on 1/4 of its obtuse faces followed by flips. The point is 
	 added at the projection of a random point arround the centroid of the face on its longest edge.
	  The first iteration of it used random points around the centroid but that was too volatile since
	 a single point inside a face creates at least 2 new obtuse faces. So after testing i decided to
	 project a single point on the longest edge. This way after the flips we usually have about the 
	 same number of obtuse faces but they are new and some of the other methods may work now.
	  After the new random steiners and flips we call a local search to try and reduce the obtuse faces.
	 Then if the obtuse faces of the copy cdt are less than the original we make them the same. This will
	 loop for a number of times until either we reached 0 obtuse faces or a certain number of loops.

Testing results:
	-Local search in most cases performed a bit worst than my previous project since it didnt allow edge
	 flips and rarely removed all the obtuse triangles.
	-Ant collony performed about the same as the local search with some cases being better and others worse
	 thanks to the randomness.
	Parameters used for testing: 
		alpha= 3.0 beta= 0.2 xi= 1.5 psi= 0.8, 0.2 kappa= 3 (if you got more threads maybe a larger number
		might be better) lambda=0.2 L=1000       
	-Simmulated annealing added  more steiner points but also reduced the obtuse angles more. I found that
	 as it looped after a point it would stop making changes to the best cdt. With an L around 1000 it would
	 most times reduce the obtuse angles to 0. The bests results (less loops needed & less steiners used) 
	 were given when used after a local search. Also the bigger the difrence between a and b (a>b) the most
	 likely it was to reduce the obtuse faces to 0. The most used method was the projection steiner and 
	 least was midpoint of longest edge which i think is expected since it is a more simplist version of
	 the projection. Also the lower the b the more time it would take since it allowed more steiner points
	 and the faces would increase making it more complex.
	 Parameters used for testing:
		alpha=3.0 beta=0.2 L=1000

