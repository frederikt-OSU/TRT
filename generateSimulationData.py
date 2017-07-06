import math

file_obj = open("simulationData.txt","w")
file_obj.write("time\t T_1\t\t T_2\t\t delta_T\n")
file_obj.close()

T_1 = [0]
T_2 = [0]
t = [0]

T_grad = 2

file_obj = open("simulationData.txt","a")
for i in range(1,100):
    
    T_1.append(290 + math.log(i) * T_grad)
    T_2.append(T_1[i] + math.sin(i * math.pi/6) * 0.19 + 4)

    file_obj.write("%d\t %f\t %f\t %f\n" % (i-1, T_1[i], T_2[i], T_2[i]-T_1[i]))

file_obj.close()


