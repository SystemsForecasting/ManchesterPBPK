# UI PBPK model in R
# Nicola Melillo, Hitesh Mistry, 07/06/2021

library(shiny)
library(shinyBS)
library(shinyjs)

shinyUI(navbarPage("Manchester PBPK",
                   
                   ### main PBPK simulation tab panel ------------------------------------------------------------------------
                   tabPanel("Inputs/Simulation",
                            sidebarLayout(
                              sidebarPanel(
                                shinyjs::useShinyjs(),
                                h4("Define model inputs for PBPK simulations"),
                                helpText("Click on the drug-specific parameter names for a brief description; click again to remove it"),
                                hr(),
                                helpText("Define schedule related parameters"),
                                fluidRow(
                                  column(6,
                                         h5('Route')
                                  ),
                                  column(6,
                                         radioButtons("Route", 
                                                      label = NULL, 
                                                      choices = list("po - solid bolus" = 2,
                                                                     "po - dissolved bolus" = 3,
                                                                     "iv - bolus" = 1,
                                                                     "iv - infusion" = 4),
                                                      selected = 2)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Infusion duration (h)')
                                  ),
                                  column(6,
                                         numericInput("inf_dur", label = NULL, value = 0.5, step=0.1, min=0, max=6)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Daily administrations')
                                  ),
                                  column(6,
                                  
                                  sliderInput("daily_admin", "",
                                              min = 1, max = 4,
                                              value = 1, step=1),
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Days of administration')
                                  ),
                                  column(6,
                                         numericInput("days", label = NULL, value = 1, min=0, max=30, step=1)
                                  )
                                ),
                                fluidRow(
                                  column(3,
                                         h5('Dose')
                                  ),
                                  column(3,
                                         radioButtons("dose_unit", 
                                                      label = NULL, 
                                                      choices = list("mg" = 1,
                                                                     "mg/kg" = 2
                                                      ),
                                                      selected = 1)
                                  ),
                                  column(6,
                                         numericInput("dose", label = NULL, value = 10)
                                  )
                                ),
                                hr(),
                                helpText("Define molecular related parameters"),
                                fluidRow(
                                  column(6,
                                         #h5('mw (g/Mol)')
                                         tipify(h5('mw (g/Mol)'),"molecular weight: it is used to calculate the diffusion coefficient in the Noyes-Whitney model describing drug dissolution in the gastro-intestinal lumen",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('MW', label = NULL, value = 500)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Type'),"type: it influences how the drug solubility in gastro-intestinal tract is calculated from the drug intrinsic solubility",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         radioButtons("type", 
                                                            label = NULL, 
                                                            choices = list("Neutral" = 0, 
                                                                           "Acid" = 1,
                                                                           "Base" = 2),
                                                            selected = 0)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('logPow'),"logarithm of the octanole to water partition coefficient: it is used to calculate the tissue to plasma partition coefficients",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('logPow', label = NULL, value = 2)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('fup'),"fraction unbound in plasma: used to calculate the tissue to plasma partition coefficients and the fraction unbound in tissues",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('fup', label = NULL, value = 0.8,min=0,max=1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         
                                         tipify(h5('BP'),"blood to plasma ratio: it is used to calculate the plasma concentration from the blood one and to calculate the tissue to plasma partition coefficient",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('BP', label = NULL, value = 0.8)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('pKa'),"acid dissociation constant: it is used to derive the solubility in the gastro-intestinal tract from the intrinsic solubility",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('pKa', label = NULL, value = 0.8)
                                  )
                                ),
                                hr(),
                                helpText("Define dissolution and absorption parameters"),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Particle Radius (um)'),"particle radius of the formulation: it is used in the Noyes-Whitney model to derive the rate constant of drug dissolution in the gastro-intestinal tract",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('r', label = NULL, value = 25)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Density of Formulation (g/L)'),"density of the formulation: it is used to calculate the hydrodynamic radius of the diffusing drug and it is used in the Noyes-Whitney model to derive the rate constant of drug dissolution in the gastro-intestinal tract",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('rho', label = NULL, value = 1000)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Intrinisic Solubility (mg/L)'),"intrinsic solubility: is the solubility of the compound in its free acid or free base form; it is used to calculate the rate constant of drug dissolution in the Noyes-Whitney model",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('Csint', label = NULL, value = 100)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Effective Permeability (10^-4 cm/s)'),"effective permeability: it is used to calculate the absorption rate constant in the small intestine compartments",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('Peff_caco2', label = NULL, value = 2)
                                  )
                                ),
                                hr(),
                                helpText("Define distribution and clearance"),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Hepatic Intrinsic Clearance (L/h)'),"hepatic intrinsic clearance: it can be defined as a proportionality constant between the drug elimination rate from the liver and the unbound concentration of drug in the liver",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('Clh', label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Renal Intrinsic Clearance (L/h)'),"renal intrinsic clearance: it can be defined as a proportionality constant between the drug elimination rate from the kidneys and the unbound concentration of drug in the kidneys",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('Clr', label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Glomerular filtration rate (GFR)'), "Glomerular filtration rate: it is the flow rate of filtered fluid through the kidneys",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         checkboxInput("GFR_flag", "add GFR", value = FALSE, width = NULL)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Enterocyte Clearance (L/h)'),"enterocyte clearance: it can be defined as a proportionality constant between the drug elimination rate from the enterocytes and the whole concentration of drug in the enterocytes",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         numericInput('Clent', label = NULL, value = 0)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         tipify(h5('Partition Coefficient Method'),"partition coefficient: for a given organ, it is defined as the ratio between the drug concentration in the tissue and in the organ outgoing blood; it is used to describe the distribution of the drugs into the various organs",placement="bottom", trigger="click")
                                  ),
                                  column(6,
                                         radioButtons("PCM", 
                                                            label = NULL, 
                                                            choices = list("Poulin & Theil" = "PT", 
                                                                           "Berezhkhovsky" = "bere"),
                                                            selected = "PT")
                                  )
                                ),
                                hr(),
                                helpText("Define specie and sex"),
                                fluidRow(
                                  column(6,
                                         h5('Species')
                                  ),
                                  column(6,
                                         radioButtons("Species", 
                                                            label = NULL, 
                                                            choices = list("Human" = "human", 
                                                                           "Mouse" = "mouse",
                                                                           "Rat" = "rat",
                                                                           "Beagle" = "beagle",
                                                                           "Dog" = "dog"),
                                                            selected = "human")
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Sex (for human only)')
                                  ),
                                  column(6,
                                         radioButtons("Sex", 
                                                      label = NULL, 
                                                      choices = list("Female" = "female", 
                                                                     "Male" = "male"),
                                                      selected = "female")
                                  )
                                ),
                                hr(),
                                helpText("Select drug specific parameters and clinical PK data and press upload"),
                                helpText("Clinical dose and route are uploaded with PK data"),
                                fluidRow(
                                  column(6,
                                         uiOutput('libraryDrugs')
                                  ),
                                  column(6,
                                         br(),
                                         actionButton("uploadDrugParam", "Upload", width='100px'),
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         uiOutput('PKData')
                                  ),
                                  column(6,
                                         br(),
                                         uiOutput('UploadPKData')
                                  )
                                ),
                                hr(),
                                helpText("Plotting options"),
                                fluidRow(
                                  column(12,
                                         checkboxInput("keepPlots", "Superimpose the simulations", value = FALSE, width = NULL)
                                  )
                                ),
                                fluidRow(
                                  column(12,
                                         checkboxInput("plotOrgansPK", "Plot PK for all the compartments", value = FALSE, width = NULL)
                                  )
                                ),
                                fluidRow(
                                  column(12,
                                         checkboxInput("logscale", "PBPK organs in semi-log10-y scale", value = FALSE, width = NULL)
                                  )
                                ),
                                helpText("Define x & y axis limits"),
                                fluidRow(
                                  column(8,
                                         uiOutput("UI_xaxis_slider")
                                  ),
                                  column(4,
                                         br(),
                                         actionButton("scalex", "scale x"),
                                  )
                                ),
                                fluidRow(
                                  column(8,
                                         uiOutput("UI_yaxis_slider")
                                  ),
                                  column(4,
                                         br(),
                                         actionButton("scaley", "scale y"),
                                  )
                                ),
                                actionButton("clearPlot", "Clear plots", width='100px'),
                                hr(),
                                helpText("Press the button to run the simulation"),
                                actionButton("runButton", "Run simulation", width='250px'),
                                hr(),
                                helpText("Perform non compartmental analysis (NCA)"),
                                helpText("For multiple administrations, NCA is performed just for the last dose"),
                                actionButton("nca", "Perform NCA", width='250px'),
                                hr(),
                                helpText("Download data & images"),
                                actionButton("download_menu", "Download menu", width='250px'),
                                hr(),
                                helpText("Press the button to reset to default values"),
                                actionButton("resetButton", "Reset", width='250px')
                              ),
                              mainPanel(
                                fluidRow(
                                ),
                                plotOutput("PK", height="1000px"),
                                br(),
                                plotOutput("PK_comp_PBPK", height="2000px"),
                                bsModal("NCAmodal", "Non Compartmental Analysis (NCA)", "nca", size = "large", dataTableOutput('table_nca')),
                                bsModal("downloadModal", "Download menu", "download_menu", size = "small", 
                                        radioButtons("radio_download_PK", label = h4("Compartments PK"), 
                                                           choices = list("only plasma concentration [mg/L]" = 1, "plasma concentration [mg/L] & PK in all the compartments [mg]" = 2),
                                                           selected = 1),
                                        downloadButton('downloadPK', 'Download PK data'),
                                        h4('Non compartmental analysis'),
                                        helpText("click 'Perform NCA' button before dowloading the results!"),
                                        downloadButton('downloadNCA', 'Download NCA'),
                                        checkboxGroupInput("checkPlots", label = h4("Plots"), 
                                                           choices = list("plasma concentration" = 1, "fraction eliminated" = 2, "fraction absorbed" = 3, "organs PK" = 4, "ACAT PK" = 5),
                                                           selected = 1),
                                        downloadButton('downloadPlots', 'Download Plots')
                                        ),
                              )
                            )
                   ),
                   ### model description tabPanel ------------------------------------------------------------------------------------
                   tabPanel("Model description",
                            shinyUI(fluidPage(
                              withMathJax(),
                              h1('Manchester PBPK Model'),
                              p(),
                              br(),
                              p('This physiologically based pharmacokinetic (PBPK) model is designed as an educational tool. The aim is to provide a simple, open source, freely downloadable PBPK model for performing basic single subject bottom-up simulations in R. The model has not been tested against a large set of compounds. All parameters were taken from the literature. The PBPK model was implemented by using the',tags$a(href="https://cran.r-project.org/web/packages/RxODE/index.html", "RxODE"),' package and non-compartmental analysis (NCA) was performed with the ',tags$a(href="https://cran.r-project.org/web/packages/PKNCA/index.html", "PKNCA"),"package."),
                              p("The model code can be found", tags$a(href="https://github.com/NicolaMelillo/ManchesterPBPK", "here"), "."),
                              hr(),
                              h3('Basic equations of the PBPK model'),
                              p('This model is composed of two parts: a PBPK model describing the distribution, metabolism and elimination of the drug in the body and a compartmental absorption and transit (CAT) based model describing events following per oral drug administration in the gut.'),
                              br(),
                              
                              h4('PBPK model for drug distribution in the body'),
                              p('The  PBPK model is an ordinary differential equations (ODE) model. This model is composed of 15 compartments, representing the lungs, brain, heart, kidneys, bones, muscles, stomach, spleen, liver, gut, pancreas, skin, fat, arterial and venous blood. The PBPK model structure is shown in the image below, were red, blue and black-dashed arrows represent arterial and venous blood flow and clearances, respectively.'),
                              div(img(src = "PBPK_model.png", height = 400), style="text-align: center;"),
                              p('The drug distributes in all the organs. Clearance is supposed to happen only in the liver and/or kidneys. The equation describing drug distribution in all the tissues except the liver, lungs, arterial and venous blood is the following.'),
                              p(withMathJax('$$\\frac{dc_i}{dt} = \\frac{Q_i}{V_i} \\bigl(c_{art} - c_{i,b,out}  \\bigr)$$')),
                              p('\\(c_{art}\\), \\(c_{i}\\) and \\(c_{i,b,out}\\) are the drug arterial blood concentration, the drug concentration in the i-th organ and the drug concentration in the i-th organ outgoing blood, respectively; \\(Q_i\\) and \\(V_i\\) are the i-th organ blood flow and volume. For solving the previous differential equation, \\(c_{i,b,out}\\) can be expressed as a function of \\(c_{i}\\). In order to do such, an instantaneous equilibration between the blood and tissue concentration in the organs is generally assumed.'),
                              p(withMathJax('$$ \\frac{c_{i}}{c_{i,b,out}} = K_{i,b} $$')),
                              p('\\(K_{i,b}\\) is the tissue to blood partition coefficient. Instead of \\(K_{i,b}\\), the tissue to plasma concentration ratio, \\(K_{i,p}\\), is more commonly used in calculations. \\(K_{i,b}\\) can be expressed as a function of \\(K_{i,p}\\): \\(K_{i,b}=K_{i,p}/BP\\), where \\(BP\\) is the equilibrium blood to plasma concentration ratio. In PBPK models, blood and plasma concentrations are commonly assumed to be in instantaneous equilibrium. \\(c_{i,b,out}\\) can now be expressed as a function of \\(c_{i}\\) and the equation of the organ distribution in the PBPK model can be expressed as follows.'),
                              p(withMathJax('$$\\frac{dc_i}{dt} = \\frac{Q_i}{V_i} \\biggl(c_{art} - \\frac{c_i}{K_{i,p}/BP}  \\biggr)$$')),
                              p('Equations for lungs (l), liver (liv), kidneys (kid), arterial blood (art) and venous blood (ven) are the following.'),
                              p(withMathJax('$$\\frac{dc_l}{dt} = \\frac{Q_{all}}{V_l} \\biggl(c_{ven} - \\frac{c_l}{K_{l,p}/BP}  \\biggr)$$')),
                              p(withMathJax('$$\\frac{dc_{liv}}{dt} = \\frac{Q_{liv,art}}{V_{liv}} \\biggl(c_{art} - \\frac{c_{liv}}{K_{liv,p}/BP}  \\biggr) + \\frac{1}{V_{liv}}\\sum_{j\\in S} \\biggl( Q_j \\frac{c_j}{K_{j,p}/BP} \\biggr) - \\frac{CL_{h,int} \\cdot fu_p}{K_{liv,p}} \\cdot \\frac{c_{liv}}{V_{liv}} + input_{GI}$$')),
                              p(withMathJax('$$\\frac{dc_{kid}}{dt} = \\frac{Q_{kid}}{V_{kid}} \\biggl(c_{art} - \\frac{c_{kid}}{K_{kid,p}/BP}  \\biggr) - \\frac{CL_{r,int} \\cdot fu_p}{K_{kid,p}} \\cdot \\frac{c_{kid}}{V_{kid}} - GFR \\cdot fu_p \\cdot \\frac{c_{art}}{V_{kid} \\cdot BP} $$')),
                              p(withMathJax('$$\\frac{dc_{art}}{dt} = \\frac{Q_{all}}{V_{art}} \\biggl( \\frac{c_l}{K_{l,p}/BP} - c_{art}  \\biggr)$$')),
                              p(withMathJax('$$\\frac{dc_{ven}}{dt} = \\frac{1}{V_{ven}}\\sum_{j\\in D} \\biggl( Q_j \\frac{c_j}{K_{j,p}/BP} \\biggr) - \\frac{Q_{all}}{V_{ven}} c_{ven} $$')),
                              p('\\(S\\) is the set of the splanchnic organs (stomach, spleen, gut, pancreas); \\(D\\) is the set of the organs whose venous output enters directly the venous blood compartment (brain, heart, kidneys, bones, muscles, liver, skin, fat); \\(fu_p\\) is the drug fraction unbound in blood; \\(Q_{all}\\) is the cardiac output; \\(Q_{liv,art}\\) is the liver arterial blood flow; \\(CL_{h,int}\\) and \\(CL_{r,int}\\) are the hepatic and renal intrinsic clearances; \\(GFR\\) is the glomerular filtration rate; \\(input_{GI}\\) is the input from the CAT based model.'),
                              br(),
                              
                              h4('Compartmental absorption and transit based model'),
                              p('Similarly to the PBPK model for drug distribution, the CAT based model is a compartmental model composed of an ODE system. The model represents the dissolution and transit out of the stomach and dissolution, transit and absorption occurring in the small intestine. In this model, the small intestine is divided into 6 different sections: one for the duodenum, two for the jejunum and three for the ileum. The CAT based model structure is shown in the image below, where St, D, J, I and LI stand for stomach, duodenum, jejunum, ileum and large intestine and \\(M_{0}\\) is the drug dose.'),
                              div(img(src = "CAT_based_model.png", height = 300), style="text-align: center;"),
                              p('Equation for describing drug transit and dissolution in the stomach are reported below.'),
                              p(withMathJax('$$ \\frac{dx_{st,s}}{dt} = -k_{t,0} x_{st,s} - k_{d,st} x_{st,s} $$')),
                              p(withMathJax('$$ \\frac{dx_{st,d}}{dt} = -k_{t,0} x_{st,d} + k_{d,st} x_{st,s} $$')),
                              p('\\(x_{st,s}\\) and \\(x_{st,d}\\) are the solid and dissolved drug amount in the stomach; \\(k_{t,0}\\) is the rate constant for drug output from the stomach and is calculated as the inverse of the gastric emptying time; \\(k_{d,st}\\) is the drug dissolution rate constant in the stomach and is described with the Noyes-Whitney model, as shown below.'),
                              p(withMathJax('$$ k_{d,st} = \\frac{3D}{\\rho h r} \\biggl( C_{st} - \\frac{x_{st,d}}{V_{st}} \\biggr) $$')),
                              p('\\(\\rho\\) is the density of the drug particle; \\(r\\) is the particle radius of the formulation; \\(h\\) is the effective thickness of the hydrodynamic diffusion layer and is calculated from \\(r\\) with the Hintz and Johnson model: \\(h=r\\) if \\(r<30 \\mu m\\), otherwise \\(h=30\\mu m\\). \\(C_{st}\\) is the drug solubility in the stomach and is calculated from the intrinsic drug solubility (\\(C_{int}\\)) as \\(C_{st}=C_{int} \\cdot \\alpha_{st} \\), where \\(\\alpha_{st}\\) is defined for neutral, monoprotic acidic and basic compounds using the Henderson Hasselbalch equation as follows.'),                              
                              p(withMathJax('$$ \\alpha_{neutral} = 1 $$')),
                              p(withMathJax('$$ \\alpha_{acid} = 1 + 10^{pH-pKa} $$')),
                              p(withMathJax('$$ \\alpha_{base} = 1 +10^{pKa-pH}$$')),
                              p('\\(pKa\\) is the drug dissociation constant while pH is the pH of the environment (either the stomach or intestine section where the drug dissolves).'),
                              p('In the Noyes-Whitney model, \\(D\\) is the drug diffusion coefficient and can be calculated from the Stokes-Einstein equation, as follows.'),
                              p(withMathJax('$$ D = \\frac{k_b T}{6\\pi\\eta_w R_h} $$')),
                              p('\\(k_b\\) is the Boltzmann constant, \\(T\\) is the absolute temperature of the body in Kelvin, \\(\\eta_w\\) is the viscosity of water at body temperature and \\(R_h\\) is the hydrodynamic radius of the diffusing drug. \\(R_h\\) is calculated as follows, assuming the drug molecule is spherical in shape.'),
                              p(withMathJax('$$ R_h = \\sqrt[3]{\\frac{3 mw}{4\\pi N_A \\rho}} $$')),
                              p('\\(mw\\) is the compound molecular weight while \\(N_A\\) is the Avogadro\'s number.'),
                              p('Equations for describing the transit, dissolution and absorption happening in the i-th section of the small intestine are reported below.'),
                              p(withMathJax('$$ \\frac{dx_{i,s}}{dt} = k_{t,i-1} x_{i-1,s} - k_{t,i} x_{i,s} - k_{d,i} x_{i,s} $$')),
                              p(withMathJax('$$ \\frac{dx_{i,d}}{dt} = k_{t,i-1} x_{i-1,d} - k_{t,i} x_{i,d} + k_{d,i} x_{i,s} - k_{a,i} x_{i,d}$$')),
                              p(withMathJax('$$ \\frac{dx_{i,ent}}{dt} = k_{a,i} x_{i,d} - CL_{ent} \\frac{x_{i,ent}}{V_{i,ent}} - Q_{i,ent} \\frac{x_{i,ent}}{V_{i,ent}}$$')),
                              p('\\(x_{i,s}\\), \\(x_{i,d}\\) and \\(x_{i,ent}\\) are the amount of solid and dissolved drug in the i-th compartment of the small interstine and the amount of drug in the i-th enterocytic compartment; \\(Q_{i,ent}\\) and \\(V_{i,ent}\\) are the blood flow and volume of the i-th section of the enterocytes; \\(k_{d,i}\\) is the drug dissolution rate constant in i-th segment of the small intestine and its definition is analogous to \\(k_{d,st}\\); \\(CL_{ent}\\) is the clearance happening in the enterocytes compartments and is supposed to be equal for all the six sections. \\(k_{t,i}\\) is the transit time constant for the i-th small intestine compartment and is calculated as \\(k_{t,i} = (SITT\\cdot l_i/l_{tot})^{-1}\\), where \\(SITT\\) is the small intestinal transit time, \\(l_i\\) is the small intestine segment length and \\(l_{tot}\\) is the total length of small intestine. \\(k_{a,i}\\) is the absorption constant of the i-th compartment of the small intestine and is calculated from the effective jejunal permeability (\\(P_{eff}\\)) as \\(k_{a,i}=2 P_{eff} /R_{i}\\), where \\(R_{i}\\) is the radius of the intestinal compartment. \\(input_{GI}\\) in the PBPK equations is defined as follows.'),                            
                              p(withMathJax('$$ input_{GI} = \\sum_{i=1}^{6} Q_{i,ent} \\frac{x_{i,ent}}{V_{i,ent}}$$')),
                              hr(),
                              
                              h3('PBPK input parameters'),
                              p('The starting point is to provide the PBPK input parameters. Currently, the units of those parameters are considered to be standard and not customizable.'),
                              tags$div(
                                tags$ul(
                                  tags$li("mw [g/mol]: molecular weight"),
                                  tags$li("Type: select if neutral, monoprotic acid or base"),
                                  tags$li("logPow: octanol to water partition coefficient, used to calculate the partition coefficients"),
                                  tags$li("fup: fraction unbound in plasma, currently used to calculate the partition coefficient"),
                                  tags$li("BP: blood to plasma ratio, used to derive the plasma concentration from blood concentration"),
                                  tags$li("pKa: dissociation constant, used to calculate water solubility in gut sections"),
                                  tags$li("r [\\(\\mu\\)m]: radius of the particle size of the formulation, used to calculate the dissolution coefficient"),
                                  tags$li("Density of the formulation [g/L]: used to calculate the dissolution coefficient"),
                                  tags$li("Intrinsic solubility [mg/L]: used to calculate the dissolution coefficient"),
                                  tags$li("Peff [\\(10^{-4}\\)cm/s]: effective permeability across the gut layer"),
                                  tags$li("Hepatic clearance [L/h]: clearance referred to total concentration of drug into the liver"),
                                  tags$li("Renal clearance [L/h]: clearance referred to total concentration of drug into the kidneys"),
                                  tags$li("Enterocyte clearance [L/h]: clearance referred to total concentration of drug into the enterocytes, assumed to be equal for all the enterocytes sections"),
                                  tags$li("Partition coefficient method: select either the Poulin & Theil or Berezhkhovsky method"),
                                  tags$li("Species: select either human, mouse, beagle or dog"),
                                )
                              ),
                              p('In addition, the user can select the desired schedule.'),
                              tags$div(
                                tags$ul(
                                  tags$li("Route: per os (po) boluses and intravenous (iv) boluses and infusions are supported; for po route, drug can be administered both in solid and dissolved form."),
                                  tags$li("Daily administrations: up to four daily doses are supported."),
                                  tags$li("Dose [mg]."),
                                  tags$li("Days: for how many days the schedule needs to be repeated.")
                                )
                              ),
                              p('All the remaining physiological parameters are fixed to mean values. Values for these parameters and relative references can be found with the source code on github, in the data directory.'),
                              hr(),
                              
                              h3('References'),
                              p("Useful references for the PBPK modelling:"),
                              tags$div(
                                tags$ul(
                                  tags$li(tags$a(href="https://doi.org/10.1038/psp.2013.41", "Jones, Rowland-Yeo 2013"), ": tutorial showing basic concepts of single subject and population PBPK models and in vitro to in vivo extrapolation."),
                                  tags$li(tags$a(href="https://doi.org/10.1002/jps.21798", "Berezhkovskiy 2009"), ": well explained theoretical bases of PBPK models."),
                                  tags$li(tags$a(href="https://doi.org/10.1016/S0378-5173(99)00147-7", "Yu Amidon 1999"), ": explanation of the compartmental absorption and transit model."),
                                  tags$li(tags$a(href="https://doi.org/10.1016/j.xphs.2018.10.033", "Grimstein et al. 2019"), ": FDA review on regulatory use of PBPK models."),
                                  tags$li(tags$a(href="https://doi.org/10.1007/s10928-018-9615-8", "Melillo et al. 2018"), ": publication from University of Manchester and University of Pavia groups considered as a reference for the equations of the CAT based model."),
                                )
                              ),
                              p("Regulatory guidances on the use of PBPK models can be found here."),
                              tags$div(
                                tags$ul(
                                  tags$li(tags$a(href="https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-reporting-physiologically-based-pharmacokinetic-pbpk-modelling-simulation_en.pdf", "European Medicines Agency (EMA), Committee for Medicinal Products for Human Use (CHMP), 2018")),
                                  tags$li(tags$a(href="https://www.fda.gov/media/101469/download", "Food and Drug Administration (FDA), Center for Drug Evaluation and Research (CDER), 2018")),
                                  tags$li(tags$a(href="https://www.oecd.org/chemicalsafety/risk-assessment/guidance-document-on-the-characterisation-validation-and-reporting-of-physiologically-based-kinetic-models-for-regulatory-purposes.pdf", "Organisation for Economic Co-operation and Development (OECD), 2021")),
                                  tags$li(tags$a(href="https://www.who.int/ipcs/methods/harmonization/areas/pbpk_models.pdf?ua=1", "International Programme on Chemical Safety (IPCS), World Health Organization (WHO), 2010")),
                                  )
                              ),
                              hr(),
                              h3('Unordered list of possible future features'),
                              tags$ul(
                                tags$li("Nonlinear metabolism in the liver."),
                                tags$li("Expansion of compounds library."),
                                tags$li("Implementation of local and global sensitivity analysis."),
                                tags$li("Permeability limited model for the liver."),
                              ),
                              hr(),
                              
                            )),
                            
                            
                            
                   ),
                   tabPanel("Authors",
                            shinyUI(fluidPage(
                              h2('Authors'),
                              tags$div(
                                tags$ul(
                                  tags$li("Nicola Melillo, University of Manchester (",tags$a(href="https://www.linkedin.com/in/nicola-melillo-4868ba107/", "LinkedIn"),"|",tags$a(href="mailto:nicola.melillo01@gmail.com", "email"),"): development of PBPK software and Shiny R app."),
                                  tags$li("Hitesh Mistry, University of Manchester (",tags$a(href="https://www.linkedin.com/in/hitesh-mistry-1ba60121/", "LinkedIn"),"|",tags$a(href="mailto:hitesh.b.mistry@gmail.com", "email"),"): conceptualization and supervision of the work, development of Shiny R app."),
                                )
                              ),
                              p('We want to thank Professor Leon Aarons of University of Manchester and Dr. Adam Darwich of KTH for the extremely valuable help in designing the app and for the critical revision.'),
                              p('The first version of this app was developed at the University of Manchester, UK.'),
                              hr(),
                              br(),
                              br(),
                            )),
                   )
)
)

