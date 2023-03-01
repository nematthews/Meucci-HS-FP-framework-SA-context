File Structure:
1. Main_Call_file.mlx executes all components of the project, excluding data processing scripts
2. Main_Project_Fns/ directory which houses reusable functions.
	- HRP_Fns
	- MV_Fns
	- HS_FP_Fns
	- Backtest_Fns3. Data/ directory to house dummy data, raw data and derived datasets.
4. Data_Retrieval/ scripts to automate data retrieval from BB for reproducible research5. A Data_Processing/ directory to house data processing scripts.• A trials/ directory houses algorithms.• A results/ directory houses the result of the applications of algorithms housed in trials/ onthe datasets housed in data/.
