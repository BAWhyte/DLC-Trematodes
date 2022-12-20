# DLC-Trematodes

<img align="right" src="https://user-images.githubusercontent.com/15988774/187104924-b8f2c25f-422c-4617-9067-1c42e431bffa.gif" width="280px">

Posture tracking of *Himasthla rhigedana* worms using DeepLabCut (DLC). Deep learning and the visualization of three worms below were done using DLC and Python.
The animated gifs on this page showing the general movement analysis of one body part done from each subject. 

The file shared in this repository is the R code for these visualizations. In summary, the code does the following:
+ Build a function for converting time data into D HH:MM:SS format
+ Transforms a data frame using pivot_longer and filters predictions with low likelihood
+ Uses ggplot stat_density_2d to create heat maps of body part positions for the entire recorded video
+ Uses gganimate to visualize the movement of one body part as simulation of the video with tracking points shown, or as a line graph showing cumulative movement

<br clear="left"/>

![HIMAS_scatter](https://user-images.githubusercontent.com/15988774/187309543-ce6dc4db-24ee-4605-b4fb-f9a3f71e456b.gif)
![HIMAS_Line2](https://user-images.githubusercontent.com/15988774/187309544-0232a5c3-4bd1-4ad6-b1cc-96f5f080ec46.gif)
