from fia import read_fia, clip_fia, area, custom_pse
import pandas as pd 

pd.set_option('display.float_format', '{:.2f}'.format)


db = read_fia("./tests/fiadb/MT")
db_current = clip_fia(db)

# 2. Call area() directly
res = area(db_current, land_type="forest", grp_by=['ADFORCD'])
print(res)
print("Rows:", len(res))

# Try our custom_pse 
"""Integration test: custom_pse on a simple forest/non-forest indicator."""
cond_list = area(
    db_current,
    land_type="forest",
    cond_list=True,
    grp_by=["ADFORCD"],
)

# TODO: The PERC_AREA_SE estimate doesn't exactly match rFIA. Fix this. 

cond_meets = cond_list.assign(
    EVAL_TYP="CURR",
    forest_flag=lambda d: (d["PROP_FOREST"] > 0).astype(float),
)

custom_pse_check = custom_pse(
    db=db_current,
    x=cond_meets,
    x_vars=["PROP_FOREST"],
    x_grp_by=["ADFORCD"],
    method="TI",
    totals=True,
    variance=True,
)

print(custom_pse_check)