# Electronic Notebook  
## Lucas Paul  

**6-23-2025**  

Working to implement changes onto pressure for better visualization 
and implement cell death aspect
Cell death based on cell-cell heavy paper via apoptosis: 2 ways to achieve
one is random cell death (highly unlikely) or critical pressure from force 
interactions between cells or fields. 
Most changes made to RadGrowth. Cell death is based on the key idea 'exempt' 
which exempts the dead cell from any interaction, treating it as dead.

**6-24-2025**  

Debugging of cell death and dynamic RGB changes in cells. Issue was in
simple changes in RadGrowth file in the access of which cells were able to change the
RGB values due to a faulty if/else statement.  
Next implementation is going to be debugging further and more rigorous tets in death/RGB pressure
as well as optimizing for computational efficiency in various files and
working to impliment better mitosis dynamics. 

**6-25-2025**

General Debugging; ready for a larger simulation to run with all other types of parameters experiencing no bugs.  
Planning to pursue accurate, large amounts of data and plots with simulation.

**6-26-2025**

Optimizing code for plotting performance. Generating and replicating plots from experimentalists
paper. Code was rewritten for usage and for optimized ease of access. Boolean conditionals were implemented
for ease of producing plots based on the visualize function. Github push not working, hoping
push will work for tomorrow. 

**6-27-2025**

Coding being used for general plotting and data analysis. Opportunities to look into next:  
time-based cell death likelihood, quantitative description of cell radius in clusters, larger box 
histograms to represent correlation between field strength and mitosis.  

 
**6-30-2025**

Implemented a few changes- worked to create a cell life length correlation to death rather than
random chance on each iteration (new cells should have a lower chance of dieing). Debugged growth
rate only being applied to a single cell rather than all cells. Allowed for better control of various 
factors of the model using logical on/off type switches. Implemented cell division based on experiment
conducted from Zhao 1999 paper where it is stated that the cells will tend to split perpendicular to the electric field, which
was implemented based on crossing cosine with y and field angle and sine with x to get the cross
splitting (perpendicular to field). 

**7-1-2025**

New repository created due to MATLAB error in pushing changes to old repository 'computational biophysics'.  
Headed towards final implementations to model and then working to replicate Abijeet's data from paper.  
Goal is to work on a decreasing velocity factor based on a function as a cell gets close to splitting.  
Both implemented, working next towards attaining more data and more accurate data for comparison.


**7-2-2025**

Working towards units of parameters and how the parameters used in this model work out to be in natural units.  
This process is completed by mathematic reasoning between the few physical constraint parameters such
as the time of simulation, cell size and active motility speed and the rest of the values can be 
determined mathematically or tweaked slightly if mathematics cannot represent it in a simple calculation.
The goal is to get the model consistant with data as well as explore new ideas while tweaking.
One idea attained was the idea of "spring force" buildup as a cell approaches mitosis, acting as a buildup of
excess repulsive force towards the two daughter cells, causing a slingshot of force away from one another.