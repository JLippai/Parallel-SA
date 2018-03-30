float channelPenalty(graph_ptr device)
{
   /* 
      As described in the Fluigi Thesis (Haiyao Huang 2016): 

      channelPenalty = channelLength*channelCost + numPenalty*overlapCost

      In other words, the total cost incurred in channel placment is the the sum of the 
      total cost of the channel lengths, plus the cost due to any overlapping channels. 

      The cost of a channel is given by the total length of channels, weighted by the cost per unit length. This is 2, as cited by Huang. 
      The cost of a overlap penalty is, reasonably, very large. Each overlap has a cost of 10000, as cited by Huang. 

      This function realizes the above. 
   */

    float cumulative_cost = 0; 
    float single_channel_cost = 0; 
    float single_overlap_cost = 0;
    float cost_per_unit_length = 2;
    int i,j,k;
    char meport;
    char youport;

    for (i = 0; i < N_COMPONENTS; i++)                              // Iterate over all components in device. 
    {
        for (j = 0; j < 4; j++)                                     // For each component, iterate over all ports. 
        {
            meport  = '';
            youport = '';
            single_channel_cost = 0;
            single_overlap_cost = 0;

            if ( (device->components[i]->neighbors[j]) > i )        // Check to make sure we have not already counted this components channel contribution.
            {
                int you_id = device->components[i]->neighbors[j];   // neighbors component id (to lookup/refference component list)
                int me_id  = i;                                     // current component id   (to lookup/reffernece component list)
                component me  = device->components[i];              // Current component 
                component you = device->components[you_id];         // Current neighbor component 

                switch (j)                                          // Find location of current components port. 
                {
                  case 0: // Top port 
                  {
                    int me_x = me.x + (me.width / 2);
                    int me_y = me.y;
                    meport = 't';
                  }
                  break;
                  case 1: // Left port
                  {
                    int me_x = me.x;
                    int me_y = me.y + (me.height / 2);
                    meport = 'l';
                  }
                  break;
                  case 2: // Right port
                  {
                    int me_x = me.x + me.width;
                    int me_y = me.y + (me.height / 2);
                    meport = 'r';
                  }
                  break;
                  case 3: // Bottom port 
                  {
                    int me_x = me.x + (me.width / 2);
                    int me_y = me.y + me.height;  
                    meport = 'b';
                  }
                  break;
                }

                for (k = 0; k < 4; k++)                       // Find location of neighbor component port. 
                {
                  if (you->neighbors[k] == me_id)
                  switch (k)
                  {
                    case 0: // Top port 
                    {
                      int you_x = you.x + (you.width / 2);
                      int you_y = you.y;
                      youport = 't'; 
                    }
                    break;
                    case 1: // Left port
                    {
                      int you_x = you.x;
                      int you_y = you.y + (you.height / 2);
                      youport = 'l';
                    }
                    break;
                    case 2: // Right port
                    {
                      int you_x = you.x + you.width;
                      int you_y = you.y + (you.height / 2);
                      youport = 'r';
                    }
                    break;
                    case 3: // Bottom port
                    {
                      int you_x = you.x + (you.width / 2);
                      int you_y = you.y + you.height;
                      youport = 'b';
                    }
                    break;
                  }
                }

                /* 
                    Recall: channelPenalty = channelLength*channelCost + numPenalty*overlapCost
                    Here, we calculate these seperate components.  
                */ 
                /* 
                    Calculate channel length cost: (channelLength*channelCost)
                    We measure the distance as a Manhattan Geometry, since channels cannot be routed diaganoly. 
                */
                single_channel_cost = (abs(you_x - me_x) + abs(you_y - me_y)) * cost_per_unit_length;

                /* 
                    Calculate overlap penalty cost: (numPenalty*overlapCost)
                    There are four types of overlaps, as described by Huang. These are checked here. 
                */
                if (meport == 't' && youport == 'b')
                  if (me_y > you_y)
                    single_overlap_cost = 10000;
                if (meport == 'b' && youport == 't')
                  if (me_y < you_y)
                    single_overlap_cost = 10000;
                if (meport == 'r' && youport == 'l')
                  if (me_x > you_x)
                    single_overlap_cost = 10000;
                if (meport == 'l' && youport == 'r')
                  if (me_x < you_x)
                    single_overlap_cost = 10000;

                cumulative_cost += (single_channel_cost + single_overlap_cost);
            }
        }
    }
    return cumulative_cost;
}