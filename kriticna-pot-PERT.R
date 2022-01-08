#Kritična pot po modelu PERT


#Podatki

#Uvozimo graf zapisan v csv datoteki kot seznam povezav
povezave <- read.csv(file = 'Grafi/graf1.csv')

u <- povezave$u
v <- povezave$v


#Generiramo porametre za beta porazdelitev
a <- floor(runif(length(u), 2, 9)) #vektor najkrajših možnih trajanj, leva krajišča intervala
b <- floor(runif(length(v), 15, 26)) #vektor najdaljših možnih trajanj, desna krajišča intervala

p <- runif(length(u), 1, 10) #vektor slučajnih vrednosti parametra p
q <- runif(length(v), 1, 10) #vektor slučajnih vrednosti parametra q
M <- ((p - 1)*b + (q - 1)*a)/(p + q - 2) #vektor najbolj verjetnih trajanj M

#Seznam povezav
opravila = data.frame(
  u,
  v)


#-------------------------------------------------------------------------------------------------------

#Funkcije


#Seznam sosedov iz seznama povezav
seznam_sosedov = function(opravila){
  m = max(opravila$v)
  sapply(1:m, function(x)opravila$v[opravila$u==x])}

#Predhodniki iz seznama povezav
seznam_predhodnikov = function(opravila){
  m = max(opravila$v)
  sapply(1:m, function(x)opravila$u[opravila$v==x])}

#Prehod po naprej
prehod_naprej = function(predhodniki, cas){
  n_tock = length(predhodniki)
  t = rep(0,n_tock)
  predniki_na_poti = c(0, rep(1,n_tock-1))
  for (i in 2:n_tock) {
    for (j in predhodniki[[i]])
      if(t[j] + cas[j] > t[i]){
        t[i] = t[j] + cas[j]
        predniki_na_poti[i] = j
      }
  }
  list(t, predniki_na_poti)
}

#Prehod po nazaj
prehod_nazaj = function(sosedi,t, cas){
  n_tock = length(sosedi)
  T = rep(t[n_tock],n_tock)
  for (i in (n_tock-1):1) {
    for (j in sosedi[[i]])
      if(T[j] - cas[i] < T[i])
        T[i] = T[j] - cas[i]
  }
  T
}   

# pot iz tabele predhodnikov
k_pot = function(predniki){
  n_tock = length(predniki)
  i = n_tock
  pot = c()
  while (i != 0){
    pot = c(i,pot)
    i = predniki[[i]]
  }
  pot
}

#metoda kritične poti
kriticna_pot = function(opravila, cas, varianca){
  n_tock = max(opravila$v)
  predhodniki = seznam_predhodnikov(opravila)
  sosedi = seznam_sosedov(opravila)
  naprej = prehod_naprej(predhodniki,cas)
  t = naprej[[1]]
  predniki = naprej[[2]]
  T = prehod_nazaj(sosedi, t, cas)
  list(t=t, T=T, criticalPath=k_pot(predniki), varianca=varianca)
}

#-------------------------------------------------------------------------------------------------------

#Ocene za pričakovana trajanja in variance


# Klasičen PERT
mi1 = (a + 4*M + b)/6
var1 = ((b - a)/6)^2


# 1. Modifikacija
# Uporabna, če za parametra, ki določata beta porazdelitev velja zveza:
# p = 1 + 2*(M - a)/(b - a) in q = 1 + 2*(b - M)/(b - a)
mi2 = (a + 2*M + b)/4
var2 = ((b + 2*M - 3*a)*(3*b - 2*M - a))/80


# 2. Modifikacija
# Uporabna, ko je znana pričakovana vrednost porazdelitve in je enaka mi1
# tedaj velja zveza: p = 1 + 4*(M - a)/(b - a) in q = 1 + 4*(b - M)/(b - a)
mi3 = mi1
var3 = ((b + 4*M - 5*a)*(5*b - 4*M - a))/252


# 3. Modifikacija
# Uporabna, ko je znana varianca porazdelitve in je enaka V, navadno vzamemo
# To nam predstavi kubično enačbo, katere edino realno rešitev označimo z qV
# Najpogosteje vzamemo V = var1
qV = runif(1, 1, 10) #nevem a tu vzamemo kar eno naključno število
mi4 = (qV*M*(b - a) + b*(b + a -2*M))/(qV*(b - a) + b + a - 2*M)
var4 = var1


# 4. Modifikacija
# Primer s poljubno porazdelitvijo z znanim M
mi5 = mi2
var5 = ((a - b)^2 + 6*(a - M)^2 + 6*(b - M)^2)/144


# 5. Modifikacija
# Primer z znano verjetnostjo na poljubnem podintervalu
# Na nekem podintervalu [a1, b1] intervala [a, b] bo verjetnost E
a1 = runif(length(u), a, b)
b1 = runif(length(v), a1, b)
E = runif(length(u), 0 , 1)  #nevem ali je vredu tako ali bi bilo potrebno generirat bolj pomenljive intervale
mi6 = E/2*(a1 + b1) + ((1 - E)/2)*(b^2 - a^2 - b1^2 + a1^2)/(b - a - b1 + a1)
var6 = ((1 - E)*(b^3 - a^3 - b1^3 + a1^3))/(3*(b - a - b1 + a1)) + E*(b1^2 + a1^2 + a1*b1) - mi6^2
#tukaj dobimo zelo velike variance


# Model CMP za primerjavo
cas = M


#-------------------------------------------------------------------------------------------------------

# Izračuni na primerih

pert = kriticna_pot(opravila, mi1, var1)
cpm = kriticna_pot(opravila, cas, rep(0, length(M)))  
pert_1m = kriticna_pot(opravila, mi2, var2)
pert_2m = kriticna_pot(opravila, mi3, var3)
pert_3m = kriticna_pot(opravila, mi4, var4)
pert_4m = kriticna_pot(opravila, mi5, var5)
pert_5m = kriticna_pot(opravila, mi6, var6)


# Kritične poti
kp_pert = pert$criticalPath
kp_cpm = cpm$criticalPath
kp_pert_1m = pert_1m$criticalPath
kp_pert_2m = pert_2m$criticalPath
kp_pert_3m = pert_3m$criticalPath
kp_pert_4m = pert_4m$criticalPath
kp_pert_5m = pert_5m$criticalPath

# Trajanja kritičnih poti
trajanje_pert = pert$t[length(pert$t)]
trajanje_cpm = cpm$t[length(cpm$t)]
trajanje_pert_1m = pert_1m$t[length(pert_1m$t)]
trajanje_pert_2m = pert_2m$t[length(pert_2m$t)]
trajanje_pert_3m = pert_3m$t[length(pert_3m$t)]
trajanje_pert_4m = pert_4m$t[length(pert_4m$t)]
trajanje_pert_5m = pert_5m$t[length(pert_5m$t)]

# Variance kritičnih poti
varianca_kriticne_poti = function(pert){
  varianca_kp = 0
  for (i in pert$criticalPath) {
    varianca_kp = varianca_kp + pert$varianca[i]
  }
  varianca_kp #varianca kriticne poti
}

var_pert = varianca_kriticne_poti(pert)
var_pert_1m = varianca_kriticne_poti(pert_1m)
var_pert_2m = varianca_kriticne_poti(pert_2m)
var_pert_3m = varianca_kriticne_poti(pert_3m)
var_pert_4m = varianca_kriticne_poti(pert_4m)
var_pert_5m = varianca_kriticne_poti(pert_5m)

