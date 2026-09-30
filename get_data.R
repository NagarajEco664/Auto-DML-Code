# get_data.R
# Treatment = log income. Price is a control. Distance driven is NOT used
# (it is a channel through which income affects gas, so controlling for it
# would give the effect holding driving fixed, not the total effect).

get_data<-function(spec,quintile){

  df  <- read.dta("gasoline_final_tf1.dta")

  N=nrow(df)

  # take logs of continuous vars
  cols.to.log<-c("gas","price","income","age")
  df[,cols.to.log]<-lapply(df[,cols.to.log],log)

  # output
  Y=df$gas

  # construct vars
  df$age2<-((df$age)^2)
  df$income2<-((df$income)^2)
  df$price2<-((df$price)^2)

  # ensure diff memory location
  df.up<-data.frame(df)
  df.down<-data.frame(df)

  # shift LOG INCOME (the treatment)
  incomes<-df$income
  delta=sd(incomes)/4          # full gap between up and down

  df.up$income<-incomes+delta/2
  df.down$income<-incomes-delta/2

  # recompute every term that depends on income
  df.up$income2<-(df.up$income)^2
  df.down$income2<-(df.down$income)^2

  # specification: income and income2 interact with the controls; price is a control
  if (spec==1){
    formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+income+income2+
      income:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+
      income2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+urban+youngsingle+
      age+age2 + price + price2
  } else {
    formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+income+income2+
      income:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+price+price2)+
      income2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+price+price2)+
      urban+youngsingle+
      age+age2 + price + price2
  }

  regressors<-model.matrix(formula,data=df)
  regressors.up<-model.matrix(formula, data=df.up)
  regressors.down<-model.matrix(formula, data=df.down)

  if (quintile>0){
    q <- ntile(df$income, 5)
    Y.q=Y[q==quintile]
    regressors.q=regressors[q==quintile,]
    regressors.up.q=regressors.up[q==quintile,]
    regressors.down.q=regressors.down[q==quintile,]
  } else {
    Y.q=Y
    regressors.q=regressors
    regressors.up.q=regressors.up
    regressors.down.q=regressors.down
  }

  return(list(Y.q,regressors.q,regressors.up.q,regressors.down.q,delta))

}
