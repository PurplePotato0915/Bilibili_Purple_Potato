search_GenAI=function(token){
  library(httr)
  library(jsonlite)
  
  url <- "https://api.siliconflow.cn/v1/chat/completions"
  messagetmp <- token
  
  system_prompt <- '
# 角色定义
role: "Assistant"
author: "DeepSeek"
'
  payload <- list(
    messages = list(
      list(role = "system", content = system_prompt),
      list(role = "user", content = messagetmp)
    ),
    model = "deepseek-ai/DeepSeek-V3"
  )
  
  headers <- add_headers(
    "Authorization" = paste("Bearer", "sk"),
    "Content-Type" = "application/json"
  )

  response <- POST(
    url,
    config = headers,
    body = toJSON(payload, auto_unbox = TRUE),  
    encode = "json"
  )
  
  content <- content(response, "parsed")
  full_content <- content$choices[[1]]$message$content
  cat(full_content)
}
