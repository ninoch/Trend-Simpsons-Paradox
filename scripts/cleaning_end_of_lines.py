import pandas as pd


# This is an example of how you can make sure your .csv file is a good .csv file!
# Please read this post: https://www.loginradius.com/engineering/eol-end-of-line-or-newline-characters/

df = pd.read_csv("../input/bad_end_of_lines.csv")
df.to_csv("../input/fine_end_of_lines.csv", line_terminator='\n', index=False)
