library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(countrycode)
library(ggthemes)
library(viridis)
library(splitstackshape)
library(ggpubr)

rm(list = ls())

terr_group<-read.csv("terrorist-groups.csv")
terr_group$EndYear<-ifelse(terr_group$EndYear==".", 2019, terr_group$EndYear)
terr_group<-terr_group%>%mutate(duration=as.numeric(EndYear)-StrYear)%>%
  mutate(terronset=1)%>%rename(year=StrYear)

duration<-ggplot(terr_group, aes(x=duration))+geom_histogram(aes(y=..density..), color="grey", fill="#E3E3E3", alpha=0.5)+
  geom_histogram(alpha=.2, color="#A0A0A0", bin = 10)+
  theme_bw()+
  labs(title="Historical Terrorist Groups Duration", x="Duration (years)", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("duration.png",  width=12.5, height=8.25, dpi=300)

terr_group$Ideology<-factor(terr_group$Ideology, levels=c("Nationalist", "Leftist", "Anarchist", "Rightist", "Religious", "Other"))
  
ideology<-ggplot(terr_group, aes(x=Ideology, fill=Ideology))+
  geom_bar(alpha=0.5)+theme_bw()+scale_fill_viridis(discrete = T)+
  labs(title="Historical Terrorist Groups Ideology", x="Ideology", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("ideology.png", width=12.5, height=8.25, dpi=300)

pre1945<-terr_group%>%filter(year<=1945)
post1945<-terr_group%>%filter(year>=1946)

ideology_pre1945<-ggplot(pre1945, aes(x=Ideology, fill=Ideology))+
  geom_bar(alpha=0.5)+theme_bw()+scale_fill_viridis(discrete = T)+
  labs(title="Historical Terrorist Groups Ideology Pre-1945", x="Ideology", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("ideology_pre1945.png", width=12.5, height=8.25, dpi=300)

ideology_post1945<-ggplot(post1945, aes(x=Ideology, fill=Ideology))+
  geom_bar(alpha=0.5)+theme_bw()+scale_fill_viridis(discrete = T)+
  labs(title="Historical Terrorist Groups Ideology Post-1945", x="Ideology", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("ideology_post1945.png", width=12.5, height=8.25, dpi=300)

onset_over_time<-ggplot(terr_group, aes(x=year, fill=terronset))+
  geom_bar(alpha=0.5)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(title="Historical Terrorist Groups Onset Over Time", x="Year", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("onset_over_time.png", width=12.5, height=4.125, dpi=300)

operation<-terr_group%>%complete(year=seq(min(1860), max(1969), by=1))
operation$duration[is.na(operation$duration)]<-0
operation<-expandRows(operation, "duration", drop = FALSE)
operation<-operation%>%rename(StrYear=year)
operation<-operation%>%
  group_by(Group)%>%
  mutate(id=row_number())%>%
  mutate(year=StrYear+id-1)
operation<-operation%>%group_by(year)%>%
  mutate(num_group=n())

operation_over_time<-ggplot(operation, aes(x=year, y=num_group))+geom_line(linetype="dashed", colour="#6E7582", alpha=0.7)+
  theme_bw()+xlim(1860, 1969)+
  labs(title="Historical Terrorist Groups in Operation Over Time", x="Year", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("operation_over_time.png", width=12.5, height=4.125, dpi=300)

theme_set(theme_pubr())

figure<-ggpubr::ggarrange(onset_over_time, operation_over_time, nrow = 2)
ggsave("combined_over_time.png", width=12.5, height=4.125, dpi=300)


##R&R Revise Figure 3##
ideology_pre1945<-ggplot(pre1945, aes(x=Ideology, fill=Ideology))+
  geom_bar(alpha=0.5)+theme_bw()+scale_fill_viridis(discrete = T)+
  labs(x=" ", y=" ")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ideology_post1945<-ggplot(post1945, aes(x=Ideology, fill=Ideology))+
  geom_bar(alpha=0.5)+theme_bw()+scale_fill_viridis(discrete = T)+
  labs(x=" ", y=" ")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

figure<-ggpubr::ggarrange(ideology_pre1945, ideology_post1945, nrow = 1, 
                          common.legend = TRUE, 
                          labels = c("Pre-1945", "Post-1945"),
                          legend = "right",
                          hjust = c(-10, -9),
                          vjust = 3,
                          font.label = list(size=8, face ="italic"))
annotate_figure(figure, top = text_grob("Historical Terrorist Group Ideology", face = "bold"),
                left = text_grob("Count"), bottom = text_grob("Ideology"))

ggsave("combined_ideology.png", width=12.5, height=4.125, dpi=300)
