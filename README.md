# Nutrient-Allocation-in-Leaf-Cutter-Ants
Designed Markovian Model to Help Interpret Experimental Results of Leafcutter Nutrient Experiment

## Table of Contents

* Supporting Documentation
  - Model Description
  - Justification for Softmax Function
* Data
  - Experimental Data
  - Simulation Data
* Scripts 
  - Simulation Script
  - Data Analysis Script

## Background

Leafcutter ants forage for leaves so that they can grow fungus gardens in their nest. These fungi require a specific ratio of carbohydrates to protein in order to grow properly. This optimal ratio is called the intake target. Any single food source will not have the ideal ratio though, so foragers must find the right combination of food items to achieve this ratio. Without a centralized command structure though, it is unclear how individual ants decide to choose some leaves over other ones. My collaborator Dr. Nathan Smith implemented an experiment where leafcutter colonies were allowed to choose two different food items which would allow them to achieve their intake target. At the end of the experiment, he found that every colony had achieved this target. He also expected to find specialists for each food type, meaning that the number of foraging trips would be correlated with the proportion of trips to one food item or the other. He instead found that there was no correlation. Ants instead seemed to be choosing randomly. How is it that random choice could result in the colony achieving its goals? I designed an individual based discrete time Markovian model to find out. 

## Methods

* Model Design
  - Colonies can store food over time, collected from two food sources that have differing levels of protein and carbohydrates. 
  - Individuals can stay in the colony, go to one food item or the other, and then return from the food item back to the colony and increase the supply of nutrients
  - The distance between the colony's current nutrient state and the intake target is the signal that ants use to forage for one food type or another
    * Ants use a popular mechanism for specialization called a response threshold to choose which food item to go to. 
    * Each ant has two response thresholds, one for each food item. The lower the threshold, the more likely an ant will choose that food item. 
    * The signal can be either positive or negative, and the sign gives which food item the ant should go to. 
    * Typical equations for the response thrshold (the Hill function) cannot handle negative signals, so I introduce a new formulation for response thresholds derived from the study of neural networks: a softmax function. My arguments for using these types of functions over softmax functions can be found here. 
  - I track individual choices over time, seeing how often virtual ants forage and the proportion of trips they make to one food item or the other. 
  
  ## Results 
  
  * Simulation results closely follow experimental results. That is, virtual ants forage as often as real ants and make the same proportion of trips to one food item:
  
  experiment
  
  * In the simulated case, there are specialists that focus on one food item or the other, but there are also a lot of cases where ants generalize and forage from both food items, which smooths out the correlation. 
  
  * The relationship between a virtual ant's threshold for one food item or the other determines where they land on this distribution. The differences in the thresholds determine how specialized an individual is on one food item or the other, while the sum determines how often they forage.  
  
  threshold
  
  * Conversely, virtual ants that choose food items randomly create a differently-shaped distribution, and thus we suspect that the real ants are not behaving randomly: 
