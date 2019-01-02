library(tools)
library(nleqslv)
library(ggplot2)
library(data.table)
#Develop IPR Curve

q_o1 <- 500
pwf1 <- 2000
q_o2 <- 685
pwf2 <- 1500


vogel_eqs <- function(z) {
  e <- numeric(2)
  f <- numeric(2)
  q_max <- z[1]
  pr <- z[2]
  
  
  f[1] <- q_max*(1-0.2*(pwf1/pr)-0.8*(pwf1/pr)^2)-q_o1
  f[2] <- q_max*(1-0.2*(pwf2/pr)-0.8*(pwf2/pr)^2)-q_o2
  f
  
}

zstart <- c(100,100)
z1 <- nleqslv(zstart,vogel_eqs)


eq = function(x){z1$x[1]*(1-0.2*(x/z1$x[2])-0.8*(x/z1$x[2])^2)}
x=pr

df1=data.frame(pwf,q_o)

for (i in 1:11) {
  
  df1[i,1]=x
  df1[i,2]=eq(x)
  x=df1[i,1]-(z1$x[2]/10)
}

#Plot IPR
ggplot()+
  geom_line(data=df1, aes(x=q_o,y=pwf, color="black"))+
  xlab('Q (bfpd)')+
  ylab('Pwf (psi)')


#ggplot(data.frame(x=c(0,z1$x[2])), aes(x=x)) + stat_function(fun=eq, geom="line") + xlab("q_o") + ylab("pwf")

#plot IPR+TPR
#open VLP50
VLP50 <-read.csv(file="VLP50.csv")
VLP75 <- read.csv(file="VLP75.csv")
VLP100 <-read.csv(file="VLP100.csv")
VLP200 <- read.csv(file="VLP200.csv")


#regresi non linear VLP1
x1<- VLP50$ql
y1<- VLP50$pwf

#regresi
plot(x1,y1,xlab="Ql", ylab = "Pwf")
model1 <-nls(y1~alpha*log(x1,exp(1))+beta, start =list(alpha=10, beta=10))
summary(model1)
nlm_fn1 <- predict(model1, newdata=VLP50$ql)
lines(VLP50$ql, nlm_fn1, col=6, lty=2)

coef1 <-coef(model1)


#regresi non linear VLP2
x2<- VLP75$qL
y2<- VLP75$pwf

plot(x2,y2,xlab="Ql", ylab = "Pwf")
model2 <-nls(y2~alpha*log(x2,exp(1))+beta, start =list(alpha=10, beta=10))
summary(model2)
nlm_fn2 <- predict(model2, newdata=VLP75$qL)
lines(VLP75$qL, nlm_fn2, col=6, lty=2)

coef2 <-coef(model2)

#regresi non linear VLP3
x3<- VLP100$qL
y3<- VLP100$pwf

plot(x3,y3,xlab="Ql", ylab = "Pwf")
model3 <-nls(y3~alpha*log(x3,exp(1))+beta, start =list(alpha=10, beta=10))
summary(model3)
nlm_fn3 <- predict(model3, newdata=VLP100$qL)
lines(VLP100$qL, nlm_fn3, col=6, lty=2)

coef3 <-coef(model3)

#regresi non linear VLP3
x4<- VLP200$qL
y4<- VLP200$pwf

plot(x4,y4,xlab="Ql", ylab = "Pwf")
model4 <-nls(y4~alpha*log(x4,exp(1))+beta, start =list(alpha=10, beta=10))
summary(model4)
nlm_fn4 <- predict(model4, newdata=VLP200$qL)
lines(VLP200$qL, nlm_fn4, col=6, lty=2)

coef4 <-coef(model4)

#Plot IPR-TPR
ggplot()+
  geom_line(data=df1, aes(x=q_o,y=pwf, color="black"))+
  geom_line(data=VLP50, aes(x=ql,y=pwf, color="yelow"))+
  geom_line(data=VLP75, aes(x=qL,y=pwf, color="red"))+
  geom_line(data=VLP100, aes(x=qL,y=pwf,color="green"))+
  geom_line(data=VLP200, aes(x=qL,y=pwf,color="blue"))+
  xlab('Q (bfpd)')+
  ylab('Pwf (psi)')


#intersection TPR-IPR1
intersection1 <- function(z) {
  f <- numeric(2)
  x1 <- z[1]
  y1 <- z[2]
  
  
  f[1] <- z1$x[1]*(1-0.2*(y1/z1$x[2])-0.8*(y1/z1$x[2])^2)-x1
  f[2] <- coef1[1]*log(x1,exp(1))+coef1[2]-y1
  f
  
}

zstart1 <- c(100,100)
z1_1 <- nleqslv(zstart1,intersection1)

x_sect1 <-z1_1$x[1]
y_sect1 <-z1_1$x[2]

#intersection TPR-IPR2
intersection2 <- function(z) {
  f <- numeric(2)
  x2 <- z[1]
  y2 <- z[2]
  
  
  f[1] <- z1$x[1]*(1-0.2*(y2/z1$x[2])-0.8*(y2/z1$x[2])^2)-x2
  f[2] <- coef2[1]*log(x2,exp(1))+coef2[2]-y2
  f
  
}

zstart2 <- c(100,100)
z1_2 <- nleqslv(zstart2,intersection2)

x_sect2 <-z1_2$x[1]
y_sect2 <-z1_2$x[2]


#intersection TPR-IPR3
intersection3 <- function(z) {
  f <- numeric(2)
  x3 <- z[1]
  y3 <- z[2]
  
  
  f[1] <- z1$x[1]*(1-0.2*(y3/z1$x[2])-0.8*(y3/z1$x[2])^2)-x3
  f[2] <- coef3[1]*log(x3,exp(1))+coef3[2]-y3
  f
  
}

zstart3 <- c(100,100)
z1_3 <- nleqslv(zstart3,intersection3)

x_sect3 <-z1_3$x[1]
y_sect3 <-z1_3$x[2]


#intersection TPR-IPR4
intersection4 <- function(z) {
  f <- numeric(2)
  x4 <- z[1]
  y4 <- z[2]
  
  
  f[1] <- z1$x[1]*(1-0.2*(y4/z1$x[2])-0.8*(y4/z1$x[2])^2)-x4
  f[2] <- coef4[1]*log(x4,exp(1))+coef4[2]-y4
  f
  
}

zstart4 <- c(100,100)
z1_4 <- nleqslv(zstart4,intersection4)

x_sect4 <-z1_4$x[1]
y_sect4 <-z1_4$x[2]

G_inj1 <- VLP50$glr[1]
G_inj2 <- VLP75$GLR[1]
G_inj3 <- VLP100$GLR[1]
G_inj4 <- VLP200$GLR[1]

#GLPC
#GLPC subset table
GLPC_table <- data.table(x=c(G_inj1, G_inj2, G_inj3,G_inj4), y=c(x_sect1,x_sect2,x_sect3,x_sect4))

#GLPC curve
ggplot()+
  geom_line(data=GLPC_table, aes(x=x, y=y))+
  xlab('Gas Lift Inject')+
  ylab('Qliq')
