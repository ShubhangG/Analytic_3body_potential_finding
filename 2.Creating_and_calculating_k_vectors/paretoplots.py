import numpy as np 
import matplotlib.pyplot as plt 
from paretochart import pareto
#34.286591
x= np.arange(3)
labels = ["Potfit", "2-body term", "3-body term"]
F=[3.2616, 0.638344652345,0.444492609527]
E=[0.0067012,0.0037734882108,0.00356012663644]
# E.append(0.0)

# for i in x:
# 	dfnf0 = np.loadtxt("eam"+repr(i)+".force", usecols=[3,5])
# 	F.append(np.sum(dfnf0[:,0])/np.sum(np.square(dfnf0[:,1])))
# 	if i!=0:
# 		dfne0 = np.loadtxt("eam"+repr(i)+".energy", usecols=[4,5])
# 		E.append(np.sum(dfne0[:,1])/np.sum(np.square(dfne0[:,0])))

# width=0.35
# plt.bar(x, F, width, label="Forces",alpha=0.4, color= 'r')
# plt.bar([p+width for p in x], E, width,label="Energy", alpha=0.4, color = 'b')
# plt.ylabel("Error")
# plt.xticks(x+width/2,labels)
# plt.title("Residues from fitting different n body term")
# plt.legend()
# plt.savefig("barplotError.png")
# plt.show()

#pareto(F)
plt.bar(x,F)
plt.title("RMS obtained: forces")
plt.ylabel("RMS Error")
plt.xticks(x,labels)
plt.savefig("PYresidBarplotA2F.png")
plt.show()
plt.close()

plt.bar(x,E,color='red')
plt.title("RMS obtained: energy")
plt.ylabel("RMS Error")
plt.xticks(x,labels)
plt.savefig("PYresidBarplotA2E.png")
plt.show()
plt.close()