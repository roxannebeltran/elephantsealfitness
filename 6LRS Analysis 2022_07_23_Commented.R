#setup----
#load packages

lapply(
  c(
    "Hmisc",
    "ggthemes",
    "tidyverse",
    "lme4",
    "car",
    "GGally",
    "viridis",
    "ggrepel",
    "survival",
    "ICC",
    "ggfortify",
    "ggpubr"
  ),
  require,
  character.only = T
)

#setwd("c:/marm/research/other/mmammals/roxanne/lrs/")
setwd("~/Publications/Submitted/ROXANNE_LRS Phenotype MS/LRS Manuscript Workflow")
D = read.csv('BeltranHernandez2022Data.csv') #load main spreadsheet; one row per seal

#Rename columns
D$Dist = D$Dist2CoastKM
D$Mass = D$MassGain
D$Date = D$DepartureDate
D$Diet = D$PropFish
D$Depth = D$MeanDepth

#average across multiple measurements for individual seals of each explanatory variable
D2 = aggregate(
  cbind(Lifespan, LRS, Dist, Mass, Date, Diet, Depth, AgeWhenTracked, Preg) ~
    FieldID,
  data = D,
  FUN = mean,
  na.action = na.pass,
  na.rm = TRUE
)
D2$Dead = 1 #column for survival analysis that indicates whether animals were dead or censored at end of study; all dead

#calculate repeatability of metrics
ICCest(x = D$FieldID, y = D$Mass)
ICCest(x = D$FieldID, y = D$Dist)
ICCest(x = D$FieldID, y = D$Date)
ICCest(x = D$FieldID, y = D$Depth)
ICCest(x = D$FieldID, y = D$Diet)

#descriptive metrics
single = D[!duplicated(D$FieldID), ]
summary(single$Lifespan)
summary(single$LRS)
summary(single$AgeWhenTracked)

#correlation between lrs and lifespan
shapiro.test(resid(lm(LRS ~ Lifespan, data = D2))) #not gaussian; fails norm test
mLRS = glm(LRS ~ Lifespan, data = D2, family = quasipoisson)
summary(mLRS) #proper model

#Annual reproductive success model----

#Models of pupping vs mass gain w/ or without duplicates averaged, etc.
#Show animals with repeated data
D$FieldID = as.factor(D$FieldID)
D_dup = D[duplicated(D$FieldID) | duplicated(D$FieldID, fromLast = TRUE), #find duplicates (both above and below)
          c("FieldID", "Mass", "Preg", "AgeWhenTracked")] #subset column names
D_dup[order(D_dup$FieldID, D_dup$Mass),] #order by FieldID, then mass
modelmass = glm(Preg ~ Mass, data = D, family = "binomial") #make glm model model 
summary(modelmass) 

#Make model predictions
predict(modelmass, newdata = data.frame(Mass = c(200, 300)), type = "response") #predict pregnancy status from arbitrary mass gain values
cbind(200:300,predict(modelmass, newdata = data.frame(Mass = 200:300), type = "response")) #arbitrary

predict(modelmass,newdata=data.frame(Mass=c(219-(0.25*sd(D$Mass)),219+(0.25*sd(D$Mass)))),type="response") #same, but with half sd

#Survival analyses----

#All metrics except diet
c4 = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Dist + Mass + Date + Depth,
data = D2)
summary(c4)

c3dist = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Dist + Mass + Date,
data = D2)
summary(c3dist)

c3depth = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Depth + Mass + Date,
data = D2)
summary(c3depth)

c2MassDate = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Mass + Depth,
data = D2)
summary(c2MassDate)

massM = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Mass,
data = D2)
summary(massM)

AIC(c4, c2MassDate, c3dist, c3depth)

#Diet analyses - no evidence for diet effects
c5pDiet = coxph(
  Surv(
    time = AgeWhenTracked,
    time2 = Lifespan,
    event = Dead
  ) ~ Dist + Mass + Date + Depth + Diet,
  data = D2
)
summary(c5pDiet) #with diet

c1Diet = coxph(Surv(
  time = AgeWhenTracked,
  time2 = Lifespan,
  event = Dead
) ~ Diet,
data = D2)
summary(c1Diet) #diet only

#Lifetime reproductive success simulation ----

nd2 = data.frame(
  Dist = mean(D2$Dist),
  Mass = 65:380,
  Date = mean(D2$Date),
  Depth = mean(D2$Depth)
) #new data to predict

#make new data frame with survival fit object and predicted mass gain values
D3 = data.frame(summary(survfit(c4, newdata = nd2))$table, Mass = 65:380)

D3 = D3 %>% select(Mass, rmean, se.rmean.) %>% 
            mutate(Response = "LRS") #subsetting columns used

#predict mean and SE of Preg~Mass model with new data
D3$preg_mean = predict(modelmass, newdata = nd2)
D3$preg_SE = predict(modelmass, newdata = nd2, se.fit = T)$se.fit
D3$LRS = (D3$rmean - 2) * plogis(D3$preg_mean) #reproductive lifespan times the probability of pregnancy each year for that mass gain value

nd = 100000 #number of draws
for (i in 1:nrow(D3)) { #for each mass value
  LRS_draws = (rnorm(n=nd, D3$rmean[i], D3$se.rmean.[i]) - 2) * #normal distribution of mean reproductive lifespan with SE, times
                plogis(rnorm(n=nd, D3$preg_mean[i], D3$preg_SE)) #probability of pregnancy each year with SE
  D3$LRS_low[i] = quantile(LRS_draws, probs = c(0.025)) #lower CI of LRS draws
  D3$LRS_high[i] = quantile(LRS_draws, probs = c(0.975))
}

#make data.frame with good format for ggplot
D4 = data.frame(
  Y_Value = c(D3$rmean - 2, D3$LRS), #reproductive lifespan and lifetime reproductive success
  Response2 = rep(c("RLS", "LRS"), each = nrow(D3)),
  Response = "LRS",
  Low95 = c(D3$rmean - 2 * D3$se.rmean. - 2, D3$LRS_low),
  Upper95 = c(D3$rmean + 2 * D3$se.rmean. - 2, D3$LRS_high),
  Mass = rep(D3$Mass),2)

#descriptive stats... look at LRS and RLS for various mass values
subset(D4, Mass == 250)
subset(D4, Mass == 219)
subset(D4, Mass == round(219 - (0.25 * sd(D$Mass))))
subset(D4, Mass == round(219 + (0.25 * sd(D$Mass))))

#FIGURE 1 (ggpairs plot of metric correlations)----

printVar = function(x, y) {
  #custom function to print results from cor.test() with appropriate labels
  vals = cor.test(x, y,
                  method = "pearson", alternative = "greater")[c("estimate", "p.value")] #one-sided
  names(vals) = c("r = ", "P = ") #specify names of variables that are output
  paste(names(vals), signif(unlist(vals), 2), collapse = "\n") #paste names and variables
}

my_fn = function(data, mapping, ...) {
  #custom function for plotting cor.test() results
  xData = eval_data_col(data, mapping$x) #takes in x for each panel
  yData = eval_data_col(data, mapping$y) #takes in y for each panel
  
  mainCor = printVar(xData, yData) #print results from cor.test()
  res = summary(lm(yData ~ xData))$coefficients[2, 4] #look at coefficients
  pcolor = "black"
  fontface = 1 #color is black, font is regular for all points until...
  if (res < .1) {
    #for significant cells
    pcolor = "blue" #color is blue
    fontface = 2 #font face is bold
  }
  
  p = ggplot(data = data, mapping = mapping) + #make plot with results from cor.test()
    annotate(
      x = 0.5,
      y = 0.5,
      label = mainCor,
      col = pcolor,
      fontface = fontface,
      geom = "text",
      size = 4
    ) +
    ylim(c(0, 1))
}

my_fn2 = function(data, mapping, ...) {   #custom function to make regression plots
  xData = eval_data_col(data, mapping$x) #takes in x for each panel
  yData = eval_data_col(data, mapping$y)#takes in y for each panel
  res = summary(lm(yData ~ xData))$coefficients[2, 4] #look at coefficients
  p = ggplot(data = data, mapping = mapping) + geom_point() +
    if (res < .1)
      geom_smooth(method = lm,
                  fill = "blue",
                  color = "blue",
                  ...) #only add lines if significant
}

ggally_mysmooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density() + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank())
}

#function for wrappin labels
wrapper <- function(x, ...){paste(strwrap(x, ...), collapse = "\n")}

#Visualize correlations among variables
p2 = ggpairs(
  D2[, c("Depth", "Diet", "Date", "Dist", "Mass")],
  switch = "both",
  columnLabels = c(
    wrapper("Foraging Depth (m)",12),
    wrapper("Fish in Diet (%)",12),
    wrapper("Departure Date (DOY)",12),
    wrapper("Distance From Coast (km)",12),
    wrapper("Mass Gain (kg)",12)
  ),
  upper = list(continuous = my_fn),
  lower = list(continuous = my_fn2),
  diag = list(continuous = ggally_mysmooth)
) + theme_few() +
  theme(
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
p2

ggsave(
  "LRS Fig1.png",
  plot = p2,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)

# FIGURE 2 (Annual survival and reproduction ~ mass gain) -------

#Figure out which seals were weighed 2+ times and had partial pregnancy success
out=aggregate(Preg~FieldID,data=D_dup,FUN=mean) #seals weighed twice or more
seals=out[out$Preg<1,]$FieldID #seals with partial pregnancy success
D$Shape=21; D$Color="black"; D$Size=1.4; D$label=NA; #set default aesthetics for most points

allcols=c("black","grey20","grey50") #colors for repeated measures seals
for(i in 1:length(seals)){ #loop through duplicated seals to set different shapes, colors, etc
  D$Shape[D$FieldID==seals[i]]=i+21 #shape for reated measures seal
  D$Color[D$FieldID==seals[i]]=allcols[i] #color for repeated measures seal
  D$Size[D$FieldID==seals[i]]=4 #size for repeated measures seal
  D$label[D$FieldID==seals[i]&D$Preg==0]=as.character(seals[i]) #label for repeated measures seal
}


massV=seq(50,450,by=25);yrv=c(5,10,13,15) #setup inputs for prediction
nd1=data.frame(Dist=mean(D2$Dist),Mass=rep(massV,length(yrv)),Date=mean(D2$Date), #dataframe for prediction
               Depth=mean(D2$Depth),AgeWhenTracked=rep(yrv,each=length(massV)),
               Lifespan=rep(yrv+1,each=length(massV)),Dead=1,Response="Survival")
nd1$Y_Value=predict(c4,newdata=nd1,type="survival",se.fit=T)$fit #predict values
col1=c(viridis(10)[c(1,4,6,9)]) #colors for age lines
nd1 = nd1%>%   mutate(label=if_else(Mass == min(Mass), as.character(paste0("Age ",AgeWhenTracked)), NA_character_)) #labels

Dtemp = D%>% dplyr::select(Dist,Mass,Date,Depth,Preg,Size,Shape,Color,label) %>%mutate(Response="Reproduction") #re-organize
Dtemp=rename(Dtemp,Y_Value=Preg) #rename column
D5=bind_rows(Dtemp,nd1) #bind rows of survival and reproduction data
D5=D5[order(D5$Size),] #reorder by size so big points are on top

Dtemp = D%>% dplyr::select(Dist,Mass,Date,Depth,Preg,Size,Shape,Color,label) %>%mutate(Response="Reproduction") #re-organize
Dtemp=rename(Dtemp,Y_Value=Preg) #rename column
D5=bind_rows(Dtemp,nd1) #bind rows of survival and reproduction data
D5=D5[order(D5$Size),] #reorder by size so big points are on top

y_lab=c("Annual Probability of Reproduction","Annual Probability of Survival") #specify y labels
names(y_lab)=unique(D5$Response)

#set up annotations
annotations <- data.frame(
  xpos = -Inf,
  ypos = Inf,
  annotateText = c("A", "B"),
  hjustvar = 0 ,
  vjustvar = 1,
  Response = names(y_lab)
) 

p3=ggplot(D5,aes(x=Mass,y=Y_Value)) + #make composite plot of survival and reproduction
  geom_smooth(data = subset(D5, Response=="Reproduction"), aes(x=Mass,y=Y_Value),
              method = glm, method.args=list(family="binomial"),se = F,colour= "grey20",size=1.5)+
  geom_label_repel(data=D5[which(D5$Shape > 21 & D5$Y_Value==0),],aes(label=label),force=0,nudge_y=c(.3,.2,.15),
                   col=D5[which(D5$Shape > 21 & D5$Y_Value==0),"Color"],bg="white",fontface="bold")+
  geom_point(data=D5[D5$Response=="Reproduction",],size=D5[D5$Response=="Reproduction","Size"],
             bg=D5[D5$Response=="Reproduction","Color"],stroke=1,col="black",
             pch=D5[D5$Response=="Reproduction","Shape"])+
  facet_grid(rows=vars(Response),scales="free",switch="both",
             labeller = labeller(Response = y_lab))+
  geom_line(data=D5[D5$Response=="Survival",],aes(group=AgeWhenTracked,col=as.factor(AgeWhenTracked)),size=1.5)+
  scale_color_manual(values=col1)+ylim(0,1)+
  geom_label_repel(data=nd1[!is.na(nd1$label),],aes(label = label),nudge_x = -30,na.rm = TRUE,size=4,col=col1,
                   bg="white",fontface="bold")+
  theme_few()+  labs(x="",y="")+theme(strip.placement = "outside")+
  theme(legend.position="none",axis.title=element_text(size=12),axis.text=element_text(size=12),
        strip.text = element_text(size = 12))+
  xlab("Mass Gain (kg)")+ geom_text_repel(
    data = annotations,
    aes(
      x = -Inf,
      y = Inf,
      hjust = hjustvar,
      vjust = vjustvar,
      label = annotateText
    ),col="black",
    segment.color = 'transparent',
    size = 5
  );p3

ggsave(
  "LRS Fig2.png",
  plot = p3,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

#FIGURE 3 (lifetime fitness metric composite) ----

p4 = ggplot(D4, aes(x = Mass, y = Y_Value, fill = Response2)) + 
  geom_ribbon(aes(ymin = Low95, ymax = Upper95, size = 1)) +
  geom_line(aes(col = Response2), size = 1) +
  scale_fill_manual(values = c(alpha("black", .2), alpha("darkred", .2)), name = "fill") +
  scale_color_manual(values = c(alpha("black", .8), alpha("darkred", .8))) +
  theme_few() +  theme(
    axis.text.y.left = element_text(color = "grey50"),
    axis.title.y.left = element_text(color = "grey50"),
    axis.text.y.right = element_text(color = "darkred"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.x = element_text(color = "black"),
    legend.position = "none",
    axis.title = element_text(size = 10)
  ) + xlab("Mass Gain (kg)") +
  scale_y_continuous( # Features of the first axis
    limits = c(0, 12),
    breaks = seq(0, 12, 2),
    name = "Lifetime Reproductive Success",
    sec.axis = sec_axis( ~ . * 1, name = "Reproductive Lifespan", breaks =
                           seq(0, 12, 2)) # Add a second axis and specify its features
  ) +  geom_segment(aes(
    x = quantile(D$Mass, .1),
    y = 10.5,
    xend = quantile(D$Mass, .1),
    yend = 10
  ),
  arrow=ggplot2::arrow(length = unit(0.15, "cm"))) +
  geom_segment(aes(
    x = quantile(D$Mass, .3),
    y = 10.5,
    xend = quantile(D$Mass, .3),
    yend = 10
  ),
  arrow=ggplot2::arrow(length = unit(0.15, "cm"))) +
  geom_segment(aes(
    x = quantile(D$Mass, .5),
    y = 10.5,
    xend = quantile(D$Mass, .5),
    yend = 10
  ),
 arrow=ggplot2::arrow(length = unit(0.15, "cm"))) +
  geom_segment(aes(
    x = quantile(D$Mass, .7),
    y = 10.5,
    xend = quantile(D$Mass, .7),
    yend = 10
  ),
  arrow=ggplot2::arrow(length = unit(0.15, "cm"))) +
  geom_segment(aes(
    x = quantile(D$Mass, .9),
    y = 10.5,
    xend = quantile(D$Mass, .9),
    yend = 10
  ),
  arrow=ggplot2::arrow(length = unit(0.15, "cm"))) +
  annotate(
    "text",
    x = quantile(D$Mass, .1) - 50,
    y = 10.9,
    label = "Mass Gain Quantiles:",
    size = 3
  ) +
  annotate(
    "text",
    x = 230,
    y = 2.5,
    label = wrapper(
      "Physiological tipping point: small reductions in mass gain lead to large reductions in lifetime reproductive success.",
      width = 33
    ),
    size = 3.6,
    hjust = 0
  ) +
  annotate(
    "text",
    x = quantile(D$Mass, .1),
    y = 10.9,
    label = "10%",
    size = 3
  ) +
  annotate(
    "text",
    x = quantile(D$Mass, .3),
    y = 10.9,
    label = "30%",
    size = 3
  ) +
  annotate(
    "text",
    x = quantile(D$Mass, .5),
    y = 10.9,
    label = "50%",
    size = 3
  ) +
  annotate(
    "text",
    x = quantile(D$Mass, .7),
    y = 10.9,
    label = "70%",
    size = 3
  ) +
  annotate(
    "text",
    x = quantile(D$Mass, .9),
    y = 10.9,
    label = "90%",
    size = 3
  )
p4

#LRS vs longevity plot
D$JLifespan = jitter(D$Lifespan - 2, 0.9) #jitter lifepsan to avoid overlapping points
D$JLRS = jitter(D$LRS, 0.9) #jitter LRS to avoid overlapping points

p1b = ggplot(D[!is.na(D$Mass), ], aes(y = JLRS, x = JLifespan)) +  #plot LRS vs lifespan
  theme_few() + geom_point(
    fill = "grey50",
    shape = 21,
    color = "black",
    size = 3,
    alpha = 0.5
  ) +
  scale_x_continuous(name = "Reproductive Lifespan",
                     limits = c(4, 20),
                     breaks = seq(4, 20, 2)) +
  scale_y_continuous(name = "Lifetime Reproductive Success",
                     limits = c(4, 20),
                     breaks = seq(4, 20, 2)) +
  annotate("text",
           x = 18,
           y = 4,
           label = "N=63 seals") +
  annotate(
    "text",
    x = 4,
    y = 16,
    label = wrapper(
      "1:1 line represents seals that reproduce every year of their reproductive lifespan (after sexual maturity).",
      width = 35
    ),
    size = 3.6,
    hjust = 0
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "grey50"),
    axis.title.y = element_text(color = "grey50"),
    axis.text.x = element_text(color = "darkred"),
    axis.title.x = element_text(color = "darkred"),
    axis.title = element_text(size = 10)
  ) +
  geom_abline(slope = 1, intercept = 0)
p1b

D6 = data.frame(variable = c(
  rep("Reproductive Lifespan", nrow(D)),
  rep("Lifetime Reproductive Success", nrow(D))
),
value = c(D$Lifespan - 2, D$LRS))

#make figure palette
figpalette <- c("Lifetime Reproductive Success" = "grey50",
                "Reproductive Lifespan" = "darkred")

#calculate mean and sd for each box
stat_box_data <- function(y) {
  return(data.frame(
    y = max(y) + 2,
    label = paste('mean =', round(mean(y), 1), '\n',
                  'sd =', round(sd(y), 1), '\n')
  ))
}

p6 = ggplot(data = D6, aes(
  x = variable,
  y = value,
  fill = as.factor(variable)
)) +
  geom_boxplot() +
  scale_fill_manual(values = figpalette) +
  xlab("") +
  scale_x_discrete(labels = c(
    paste("Lifetime", '\n', "Reproductive", '\n', "Success"),
    paste("Reproductive", '\n', "Lifespan")
  )) +
  scale_y_continuous(name = "",
                     limits = c(4, 20),
                     breaks = seq(4, 20, 2)) +
  theme_few() + theme(
    legend.position = "none",
    strip.placement = "outside",
    axis.text.x = element_text(colour = figpalette),
    axis.text.y = element_text(colour = "black")
  ) +
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.9,
    color = figpalette,
    size = 3
  )
p6

p1d = ggarrange(
  p4,
  labels = "A",
  hjust = -4.5,
  vjust = 2.5,
  nrow = 2,
  ggarrange(
    plotlist = list(p1b, p6),
    ncol = 2,
    heights = c(1, 1),
    widths = c(1, .6),
    align = "h",
    labels = c("B", "C"),
    hjust = -4.5,
    vjust = 2.5
  )
)

ggsave(
  "LRS Fig3.png",
  plot = p1d,
  width = 7,
  height = 7,
  units = "in"
)

