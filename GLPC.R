library(tools)
library(nleqslv)
library(ggplot2)
library(data.table)
#Develop IPR Curve

q_o1 <- 500
pwf1 <- 2000
q_o2 <- 685
pwf2 <- 1500
pwf2


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

#declared the possibility of model
#model_optional
#asumsi alpha -> saat ql ~ 0
alpha <- 1800
model1_exponential <- nls(y1 ~ alpha*exp(beta*x1), start = list(alpha=alpha, beta=2*log(2)/alpha))
model1_logarithmic <- nls (y1~alpha*log(x1,exp(1))+beta, start =list(alpha=alpha, beta=2*log(2)/alpha))
model1_polynomial1 <- lm(y1~ poly(x1,2))
model1_power <- nls(y1~alpha*x1^beta, start=list(alpha=alpha, beta=2*log(2)/alpha))

#R squared calculation
R_exponential <- 1-(deviance(model1_polynomial1)/sum((y1-mean(y1))^2))
R_logarithmic <- 1-(deviance(model1_logarithmic)/sum((y1-mean(y1))^2))
R_polynomial1 <- 1-(deviance(model1_polynomial1)/sum((y1-mean(y1))^2))
R_power <- 1-(deviance(model1_power)/sum((y1-mean(y1))^2))


data_R <-c("R_exponential", "R_logarithmic", "R_polynomial1", "R_power")
value_R <-c(R_exponential, R_logarithmic, R_polynomial1,R_power)
df_R <-data.frame(data_R, value_R)

best_model <-apply(df_R, 2, function(x) max(x, na.rm = TRUE))
a<-as.character(best_model[1])
a

#build function for regression with highest R square
best_regression <- function(x) {
  if (x == "R_exponential"){
    nlm_fn1 <-predict(model1_exponential, newdata = VLP50$ql)
    lines(VLP50$ql, nlm_fn1, col=6, lty=2)
    coef1 <-coef(model1_exponential)
    
    #intersection
    intersection1 <- function(z) {
      f <- numeric(2)
      x1 <- z[1]
      y1 <- z[2]
      
      
      f[1] <- z1$x[1]*(1-0.2*(y1/z1$x[2])-0.8*(y1/z1$x[2])^2)-x1
      f[2] <- coef1[1]*exp(coef1[2]*z[1])-z[2]
      f
      
    }
    
    zstart1 <- c(100,100)
    z1_1 <- nleqslv(zstart1,intersection1)
    
    x_sect1 <-z1_1$x[1]
    y_sect1 <-z1_1$x[2]
    
    
  }else if (x =="R_logarithmic"){
    nlm_fn1 <-predict(model1_logarithmic, newdata = VLP50$ql)
    lines(VLP50$ql, nlm_fn1, col=6, lty=2)
    coef1 <-coef(model1_logarithmic)
    #intersection
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
    
  }else if (x=="R_polynomial1"){
    nlm_fn1 <-predict(model1_polynomial1)
    lines(VLP50$ql, nlm_fn1, col=6, lty=2)
    coef1 <-coef(model1_polynomial1)
    
    #intersection
    intersection1 <- function(z) {
      f <- numeric(2)
      x1 <- z[1]
      y1 <- z[2]
      
      
      f[1] <- z1$x[1]*(1-0.2*(y1/z1$x[2])-0.8*(y1/z1$x[2])^2)-x1
      f[2] <- coef1[3]*x1^2+coef1[2]*x1+coef1[1]-y1
      f
      
    }
    
    zstart1 <- c(100,100)
    z1_1 <- nleqslv(zstart1,intersection1)
    
    x_sect1 <-z1_1$x[1]
    y_sect1 <-z1_1$x[2]
    
    
  }else {
    nlm_fn1 <-predict(model1_power, newdata = VLP50$ql)
    lines(VLP50$ql, nlm_fn1, col=6, lty=2)
    coef1 <-coef(model1_power)
    
    #intersection
    intersection1 <- function(z) {
      f <- numeric(2)
      x1 <- z[1]
      y1 <- z[2]
      
      
      f[1] <- z1$x[1]*(1-0.2*(y1/z1$x[2])-0.8*(y1/z1$x[2])^2)-x1
      f[2] <- coef1[1]*x1^coef[2]-y1
      f
      
    }
    
    zstart1 <- c(100,100)
    z1_1 <- nleqslv(zstart1,intersection1)
    
    x_sect1 <-z1_1$x[1]
    y_sect1 <-z1_1$x[2]
    
  }
  
}

#input the best regression model type to function of regression
best_regression(a) 




#Plot IPR-TPR
ggplot()+
  geom_line(data=df1, aes(x=q_o,y=pwf, color="black"))+
  geom_line(data=VLP50, aes(x=ql,y=pwf, color="yelow"))+
  geom_line(data=VLP75, aes(x=qL,y=pwf, color="red"))+
  geom_line(data=VLP100, aes(x=qL,y=pwf,color="green"))+
  geom_line(data=VLP200, aes(x=qL,y=pwf,color="blue"))+
  xlab('Q (bfpd)')+
  ylab('Pwf (psi)')




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
