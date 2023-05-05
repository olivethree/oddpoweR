library(shiny)
library(lmerTest)
library(simr)
library(magrittr)
library(tidyverse)
library(vroom)
library(cowplot)
library(ggplot2)

ui <- fluidPage(
  
  # Theme
  theme = bslib::bs_theme(bootswatch = "yeti"),
  
  # App title
  titlePanel("oddpoweR: Power analysis for Generalized Linear Mixed-Effect Models", 
             windowTitle = "oddpoweR"),
  
  br(),
  markdown("Version 0.1"),
  p("Note: The current version only computes power for two or three way interactions."),
  br(),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("select_input", "Select a model:", 
                  choices = c("2w * 2w + (1|subject) + (1|trial)",
                              "2w * 2w * 2w + (1|subject) + (1|trial)")),
      
      selectInput("modelFam", "Select model family:",
                  choices = c("binomial", "poisson"),
                  selected = "binomial"),
      
      numericInput("min_samplesize", "Minimum sample size (subjects) to simulate for:", value = 30),
      numericInput("max_samplesize", "Maximum sample size (subjects) to simulate for:", value = 35),
      numericInput("step", "Step between sample sizes in power curve", value = 1),
      numericInput("ntrials", "Nr. trials per cell (e.g. how many trials in A1_B1, or A1_B1_C2, etc.)", value = 100),
      numericInput("nsims", "Nr. iterations", value = 5),
      numericInput("input_seed", "Seed number", value = 42),
      
      conditionalPanel(
        condition = "input.select_input == '2w * 2w + (1|subject) + (1|trial)'",
        numericInput("OR_A", "Odds Ratio Predictor A:", value = 2),
        numericInput("OR_B", "Odds Ratio Predictor B:", value = 2),
        numericInput("OR_AB", "Odds Ratio Interaction A:B:", value = 2),
        numericInput("mod_int", "Model intercept", value = 1),
        numericInput("rand_id", "Random intercept for subject", value = 0.5),
        numericInput("rand_trial", "Random intercept for trial", value = 0.5)
      ),
      
      conditionalPanel(
        condition = "input.select_input == '2w * 2w * 2w + (1|subject) + (1|trial)'",
        numericInput("OR_A", "Odds Ratio Predictor A:", value = 2),
        numericInput("OR_B", "Odds Ratio Predictor B:", value = 2),
        numericInput("OR_C", "Odds Ratio Predictor C:", value = 2),
        numericInput("OR_AB", "Odds Ratio Interaction A:B:", value = 2),
        numericInput("OR_BC", "Odds Ratio Interaction B:C:", value = 2),
        numericInput("OR_AC", "Odds Ratio Interaction A:C:", value = 2),
        numericInput("OR_ABC", "Odds Ratio Interaction A:B:C", value = 2),
        numericInput("mod_int", "Model intercept", value = 1),
        numericInput("rand_id", "Random intercept for subject", value = 0.5),
        numericInput("rand_trial", "Random intercept for trial", value = 0.5)
      ),
      
      actionButton("calculate", "Determine Sample Size")
      
     
    ),
    mainPanel(
      fluidRow(
        markdown("**Effect size conversion**"),
        markdown(
        " | Odds Ratio | Cohen's *d* | Size |
        |:------|:------|:--------|
        | 1.68 | 0.20 | small  |
        | 3.47 | 0.50 | medium |
        | 6.71 | 0.80 | large  |
        ")
      ),
      br(),
      verbatimTextOutput("sim_timestamp"),
      verbatimTextOutput("ss_result"),
      verbatimTextOutput("result"),
      br(),
      plotOutput("powercurve_plot"),
      br(),
      markdown("**Save results?**"),
      downloadButton('downloadData', 'Download timestamped results (.csv)'),
      br(),
      br(),
      markdown("**Credits**"),
      markdown("This app has a strong dependency on the 'simr' R package (Green & MacLeod, 2016) and is meant to facilitate its use by a wider audience."),
      br(),
      markdown("**Citation**"),
      tags$div(HTML("<p>Oliveira, M. (2023). oddpoweR: Power analysis for generalized linear mixed-effects models. R Shiny application (Version 0.1), <a href='http://olivethree.shinyapps.io/oddpoweR'>http://olivethree.shinyapps.io/oddpoweR</a></p>")),
      markdown("**References**"),
      markdown("Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio? Interpreting the magnitudes of odds ratios in epidemiological studies. *Communications in Statistics — Simulation and Computation*, *39*(4), 860-864."),
      markdown("Green, P., & MacLeod, C. J. (2016). SIMR: An R package for power analysis of generalized linear mixed models by simulation. *Methods in Ecology and Evolution*, *7*(4), 493-498.")
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$calculate, {
    
    # Configure progress bar
    progress <- Progress$new(session, min = 0, max = 1)
    on.exit(progress$close())
    progress$set(message = 'Running simulation',
                 detail = 'This may take a while...')
    
    # Initialize simulation results data frame
    sim_results <- data.frame()
    
    if (input$select_input == "2w * 2w + (1|subject) + (1|trial)") {
      
      # Minimum sample size to simulate for
      min_ss <- isolate(input$min_samplesize)
      # Maximum sample size to simulate for
      max_ss <- isolate(input$max_samplesize)
      # Step between sample sizes
      sim_step <- isolate(input$step)
      
      # Model family
      model_fam <- input$modelFam
      
      # Run simulation
      for (i in seq(min_ss, max_ss, sim_step)) {
        
        progress$inc(amount = 0.1)
        
        # Generate data
        nr_subjects <- i
        nr_trials <- isolate(input$ntrials)
        
        pred_A_levels <- c("A1", "A2")
        pred_B_levels <- c("B1", "B2")
        
        max_trial_ids <- nr_trials * (length(pred_A_levels)*length(pred_B_levels))
        trial_levels <- sapply(1:max_trial_ids, function(x) paste0("trial_",x))
        trials <- rep(1:nr_trials, times = length(pred_A_levels)*length(pred_B_levels))
        subjects <- rep(1:nr_subjects, each = max_trial_ids)
        predictor_A <- rep(pred_A_levels, each = max_trial_ids / 2)
        predictor_B <- rep(pred_B_levels, each = max_trial_ids / 4)
        
        # Data structure without DV
        covars <- data.frame(id = subjects, 
                             trial_per_condition = rep(trials, times = nr_subjects),
                             trials = rep(trial_levels, times = nr_subjects),
                             pred_A = rep(predictor_A, times = nr_subjects),
                             pred_B = rep(predictor_B, times = nr_subjects))
        
        # Slopes for predictors
        coef_A <- isolate(input$OR_A)
        coef_B <- isolate(input$OR_B)
        coef_AB <- isolate(input$OR_AB)
        
        # Model intercept
        intercept_model <- isolate(input$mod_int)
        
        # Intercept and fixed effect coefficients derived from odds ratios
        fixed <- c(intercept_model,
                   log(coef_A),
                   log(coef_B),
                   log(coef_AB))
        
        # Random intercepts for participants clustered by class
        rand_id_value <- isolate(input$rand_id)
        rand_trial_value <- isolate(input$rand_trial)
        rand <- list(rand_id_value, rand_trial_value)
        
        if(model_fam == "binomial") {
          
          # Generate data with model parameters
          fitted.model <-
            simr::makeGlmer(
              y ~ pred_A * pred_B + (1 | trials) + (1 | id),
              family = binomial(link = "logit"),
              fixef = fixed,
              VarCorr = rand,
              data = covars
            )
        }
        
        if(model_fam == "poisson") {
          
          # Generate data with model parameters
          fitted.model <-
            simr::makeGlmer(
              y ~ pred_A * pred_B + (1 | trials) + (1 | id),
              family = poisson,
              fixef = fixed,
              VarCorr = rand,
              data = covars
            )
        }
        
        
        # Run simulation iteration
        run_results <- powerSim(
          fitted.model,
          nsim = input$nsims,
          seed = input$input_seed,
          test = fcompare(y ~ pred_A * pred_B)
        )
        
        # Extract results from output and add sample size information
        out_summary <-
          summary(run_results) %>% mutate(sample_size = nr_subjects)
        
        # Attach results to final data frame
        sim_results <- sim_results %>% rbind(out_summary)
        
        # Print object
        sim_results
      }
      
      # Add more info to results
      sim_results <- sim_results %>%
        mutate(nr_stimuli = input$ntrials,
               timestamp = Sys.time(),
               power_for = "two_way_interaction",
               design = "2w * 2w + (1|subject) + (1|trials)")
      
    }
    
    if(input$select_input == "2w * 2w * 2w + (1|subject) + (1|trial)") {
      
      # Minimum sample size to simulate for
      min_ss <- isolate(input$min_samplesize)
      # Maximum sample size to simulate for
      max_ss <- isolate(input$max_samplesize)
      # Step between sample sizes
      sim_step <- isolate(input$step)
      
      # Run simulation
      for(i in seq(min_ss, max_ss, sim_step)){
        
        # progress$set(value = i)
        progress$inc(amount = 0.1)
        
        # Generate data
        nr_subjects <- i
        nr_trials <- isolate(input$ntrials)
        
        pred_A_levels <- c("A1", "A2")
        pred_B_levels <- c("B1", "B2")
        pred_C_levels <- c("C1", "C2")
        
        max_trial_ids <- nr_trials * (length(pred_A_levels)*length(pred_B_levels)*length(pred_C_levels))
        trial_levels <- sapply(1:max_trial_ids, function(x) paste0("trial_",x))
        trials <- rep(1:nr_trials, times = length(pred_A_levels)*length(pred_B_levels))
        subjects <- rep(1:nr_subjects, each = max_trial_ids)
        predictor_A <- rep(pred_A_levels, each = max_trial_ids / 2)
        predictor_B <- rep(pred_B_levels, each = max_trial_ids / 4)
        predictor_C <- rep(pred_C_levels, each = max_trial_ids / 8)
        
        # Data structure without DV
        covars <- data.frame(id = subjects, 
                             trial_per_condition = rep(trials, times = nr_subjects),
                             trials = rep(trial_levels, times = nr_subjects),
                             pred_A = rep(predictor_A, times = nr_subjects),
                             pred_B = rep(predictor_B, times = nr_subjects),
                             pred_C = rep(predictor_C, times = nr_subjects))
        
        
        # Slopes for predictors
        coef_A <- isolate(input$OR_A)
        coef_B <- isolate(input$OR_B)
        coef_C <- isolate(input$OR_C)
        coef_AB <- isolate(input$OR_AB)
        coef_AC <- isolate(input$OR_AC)
        coef_BC <- isolate(input$OR_BC)
        coef_ABC <- isolate(input$OR_ABC)
        
        # Model intercept
        intercept_model <- isolate(input$mod_int)
        
        # Intercept and fixed effect coefficients derived from odds ratios
        fixed <- c(intercept_model,
                   log(coef_A),
                   log(coef_B),
                   log(coef_C),
                   log(coef_AB),
                   log(coef_AC),
                   log(coef_BC),
                   log(coef_ABC))
        
        # Random intercepts for participants clustered by class
        rand_id_value <- isolate(input$rand_id)
        rand_trial_value <- isolate(input$rand_trial)
        rand <- list(rand_id_value, rand_trial_value)
        
        # Generate data with model parameters
        fitted.model <- 
          simr::makeGlmer(y ~ pred_A * pred_B * pred_C + (1|trials) + (1|id),
                          family=binomial(link = "logit"), 
                          fixef=fixed, 
                          VarCorr=rand,
                          data=covars)
        
        # Run simulation iteration
        run_results <- powerSim(fitted.model,
                                nsim = input$nsims,
                                seed = input$input_seed,
                                test = fcompare(y ~ pred_A * pred_B))
        
        # Extract results from output and add sample size information
        out_summary <- summary(run_results) %>% mutate(sample_size = nr_subjects)
        
        # Attach results to final data frame
        sim_results <- sim_results %>% rbind(out_summary)
        
        # Print object
        sim_results
        
      } 
      
      # Add more info to results
      sim_results <- sim_results %>%
        mutate(nr_stimuli = input$ntrials,
               timestamp = Sys.time(),
               power_for = "three_way_interaction",
               design = "2w * 2w * 2w + (1|subject) + (1|trials)")
    }
    

    # pcurve_results <- powerCurve(fitted.model, 
    #                              test=fcompare(y ~ pred_A + pred_B),
    #                              along="id",
    #                              nsim = 5,
    #                              seed = 42)
    # 
    # sim_results <- summary(pcurve_results) %>% mutate(sample_size = nr_subjects)
    
    # Download results
    if(input$select_input == "2w * 2w + (1|subject) + (1|trial)") {
      output$downloadData <- downloadHandler(
        filename = function() { 
          paste("power_results_2w2w_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(sim_results, file)
        })
    }
    
    if(input$select_input == "2w * 2w * 2w + (1|subject) + (1|trial)") {
      output$downloadData <- downloadHandler(
        filename = function() { 
          paste("power_results_2w2w2w_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(sim_results, file)
        })
    }
    
    # Display results table
    output$result <- renderPrint({
      sim_results %>% 
        dplyr::select(`Sample Size` = sample_size, Power = mean, 
                      `Lower 95% CI` = lower, `Upper 95% CI` = upper)
    })
    
    # Display time stamp 
    output$sim_timestamp <- renderPrint({
      paste("Timestamp:", Sys.time())
    })
    
    # Display sample size result
    output$ss_result <- renderPrint({
      min_ss <- sim_results %>%
        filter(lower >= 0.80) %>%
        slice(which.min(sample_size)) %>%
        pull(sample_size)
      
      paste("Estimated required sample size:", min_ss)
    })
    
    # plot power curve
    output$powercurve_plot <- renderPlot({
      
      # # Extract computed power closest to the 80% power threshold
      pwr_thresh <- sim_results %>%
        filter(mean >= 0.80) %>%
        slice(which.min(mean)) %>%
        pull(mean)
      
      # Extract sample size closest to 80% power threshold
      min_ss <- sim_results %>%
        filter(mean >= 0.80) %>%
        slice(which.min(sample_size)) %>%
        pull(sample_size)
      
      # Plot power curve
      sim_results %>%
        ggplot(aes(x = sample_size, y = mean)) +
        geom_smooth(method=lm, formula = y ~ poly(x, 3), se = FALSE, color = "goldenrod2") +
        geom_point(color = "#CD6090") +
        geom_errorbar(aes(ymin=lower, ymax=upper), width = .1) +
        geom_hline(yintercept = .80) +
        geom_hline(yintercept = pwr_thresh,  color="gray45", alpha = .5, linetype = "dotted") +
        geom_vline(xintercept = min_ss,  color="gray45", alpha = .5, linetype = "dotted") +
        labs(y = "Power (+/- 95% CI)",
             x = "Sample size",
             title = paste0("Estimated minimum required sample size: ", min_ss),
             caption = "Sample size based on power estimate for which lower bound 95% CI > 80% threshold.") +
        theme_cowplot()

    })
  })
}

shinyApp(ui = ui, server = server)