Created by: ΑΓΓΕΛΟΣ ΒΕΝΕΤΗΣ ΦΑΝΟΥΡΑΚΗΣ sdi1300022

Files and Directories:

├── build 	<-- Inside will be built the executable Triangulation
├── Header Files 	<-- Contrains the header files
│   ├── boost_utils.hpp		<-- Functions used to read from json file	
│   ├── CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2.h	<--Custom cdt class used.Added more methods to avoid flips.
│   └── triangulation_utils.hpp		<--Functions rest of functions used between the files
├── make		<-- Inside there is the CMakeLists.txt
│   └── CMakeLists.txt
└── Sources		<-- All the source code for the functions used.
    ├── boost_utils.cpp		<--Source code of functions used to read the json file.
    ├── LocalSearch.cpp		<--Source code of the local search.
    ├── main.cpp		<--Source code of the main function.
    ├── outputfile.cpp		<--Source code of the function used to create the output.
    ├── Previous_Project.cpp	<--Source code of the previous project. Used when delaunay=false
    │			   	   to recreate the triangulation produced by the previous project.
    ├── SimulatedAnnealing.cpp	<--Source code of the sumulated annealing.
    ├── SteinerPoints.cpp	<--Source code of the functions used to find and add steiner points.
    ├── triangulation_utils.cpp <--Source code of other helping functions.
    └── AntCollony		<--Source code of the ant collony method.
How to compile and run:
	Inside the make directory there is the CMakeLists.txt. If you run cmake. inside and after
	make the files will be compiled and in the build directory the executable Triangulation
	will be created.
	To run the program you call ./Triangulation -i /path/to/input.json -o /path/to/output.json

How LocalSearch works:
	Local search goes through all the obtuse faces of triangulation and for each test each of the 5
	methods to see how many obtuse faces would be after adding a steiner.
	The 5 methods are:
		-steiner on the midpoint of the longest edge
		-steiner on the projection of the obtuse angle
		-steiner on the centroid
		-steiner on the circumcenter:
			To check the steiner on the circumcenter first the circumcenter needs to be inside
			a neighbor. If it is and the face isnt constrained we remove it. Then we add the
			steiner and the vertex of the old obtuse angle as a constrain. After we add back
			the rest of the removed vertices and make constrains the edges on the boundary of
			the polygon that the original face and neighbor had. Then we remove the the constrains
			and count the obtuse faces. This way we make a local re triangulation.
		-steiner on the centroid of merged obtuse neighbors:
			For this steiner to be added we the face to have an obtuse neighbor, both of them to be
			non constrained and the merged polygon to convex. If those conditions are met we start
			by finding the centroid of the polygon. Then we remove all the vertices. Then we add 
			the steiner point on the centroid.We then add back all the vertices and make the edges
			we want to force a constrains. After that we remove the constrains and count the obtuse
			faces. We do that check for each neighbor of the face and return the one we got the best
			result from along with its index.

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
	  Also there is the improveTriangulationthat is called for each ant. improveTriangulation() chooses one
	of the obtuse faces and a type of steiner at random. After it simulates adding the steiner to the 
	triangulation and checks saves the energy that steiner would produce. After the face that was chosen
	is returned and the energy and method that was used are stored in vectors given as arguments.
	  The main function is the ant_colony. It start by creating some vectors to store the pheromones,
	ant_energy, ant_method and calculating the starting energy. Then it starts looping for L. in each
	loop a vector ant_faces is created to store the faces that each ant chose to edit. After there is
	another loop for each ant. In that loop improvedTriangulation() is called to simulate what edits each
	ant would do. The method, face and energy each ant chose are stored in the vectors. After that we
	choose the ant that improved the energy the best (if it exists) its steiner is added to the cdt and
	the pheromon of the method is updated. Then we keep choosing the next best (if it exists) and check
	if the face it wants to edit still exists if it does the same thing happens. This way we can filter
	conflicts by applying the best first.
Testing results:
	-Local search in most cases performed a bit worst than my previous project since it didnt allow edge
	 flips and rarely removed all the obtuse triangles.
	-Simmulated annealing added  more steiner points but also reduced the obtuse angles more. I found that
	 as it looped after a point it would stop making changes to the best cdt. With an L around 1000 it would
	 most times reduce the obtuse angles to 0. The bests results (less loops needed & less steiners used) 
	 were given when used after a local search. Also the bigger the difrence between a and b (a>b) the most
	 likely it was to reduce the obtuse faces to 0. The most used method was the projection steiner and 
	 least was midpoint of longest edge which i think is expected since it is a more simplist version of
	 the projection. Also the lower the b the more time it would take since it allowed more steiner points
	 and the faces would increase making it more complex.


