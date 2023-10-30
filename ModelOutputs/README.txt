InfectionTimes - outputs the times of all infection events

col 1 - time event occurs  (years)
col 2 - cell number
col 3 - x coordinate
col 4 - y coordinate
col 5 - event type 1 (1=residential, 2=commercial)
col 6 - event type 2 (1=S->E, 2=E->C, 3=C->I)


CellInfections - outputs information about all the cells and their final amount of infected citrus

col 1 - cell number
col 2 - x coordinate
col 3 - y coordinate
col 4 - amount of residential citrus (in discrete units where 1000=full density)
col 5 - amount of commercial citrus (in discrete units where 100=full density)
col 6 - time of first infection in residential citrus 
col 7 - time of first infection in commercial citrus 
col 8-10 - final amount of residential citrus in E,C,I components
col 11-13 - final amount of commercial citrus in E,C,I components