// Custom JavaScript for the application

$(document).ready(function() {
  // Enable tooltips
  $('[data-toggle="tooltip"]').tooltip();
  
  // Custom handling for follow-up questions
  $(document).on('click', '#followup1, #followup2, #followup3', function() {
    var question = $(this).text();
    $('#ai_prompt').val(question);
    $('#submit_prompt').click();
  });
  
  // Custom handling for tab changes
  $('a[data-toggle="tab"]').on('shown.bs.tab', function(e) {
    // Trigger resize events for plots when tabs are shown
    $(window).trigger('resize');
  });
  
  // Custom handling for download buttons
  $('.download-button').click(function() {
    // Add any custom download handling here
    console.log('Download initiated for: ' + $(this).attr('id'));
  });
  
  // Custom handling for the AI analysis button
  $('#run_ai_analysis').click(function() {
    // Show loading spinner
    $('#ai_summary_report').html('<div class="loading-spinner"><i class="fa fa-spinner fa-spin"></i> Generating insights...</div>');
  });
});