import os

old = "RKO"
new = "GOU"

def rename(path,old,new):
    for f in os.listdir(path):
        os.rename(os.path.join(path, f), 
                  os.path.join(path, f.replace(old, new)))


rename(".", old=old, new=new)