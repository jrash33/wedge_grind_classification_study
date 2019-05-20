# wedge_grind_classification_study

The purpose of this study is to predict optimal wedge grind for players without having to actually fit them

- This theory suggests that there is a certain swing characteristic that influences good/bad grind performance

Steps:
1. Get test data on 20-30 players
- This is to get enough swings to be able to classify 1-5 swing types (cluster players)

2. Build algorithm
- Run R code to create full swing classification variables (heirarchical clustering)
- Using K-Means clustering, transform continuous variables into cluster of 1-5 groups (TMAN/MoCap variables)
- Random Forest method chosen as most accurate method to predict "dummy" grinds

3. Fit same players to optimal grind performance
- Evaluate algorithm accuracy

